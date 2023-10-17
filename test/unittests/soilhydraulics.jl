using Test
import Emerald.EmeraldPhysics.Constant as CS
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.PhysicalChemistry as PC
import Emerald.EmeraldLand.SoilHydraulics as SH
import Emerald.EmeraldLand.SPAC

@testset verbose = true "SoilHydraulics.jl" begin
    @testset "Vulnerability curves (k, ψ, θ)" begin
        for vc in [NS.BrooksCorey{Float64}(1), NS.VanGenuchten{Float64}("Loam")]
            θs = collect(Float64, vc.Θ_RES+eps():0.01:vc.Θ_SAT+0.02);
            ps = SH.soil_ψ_25.((vc,), θs);
            @test all(ps .<= 0);

            ps = 0:-0.1:-10;
            θs = SH.soil_θ.((vc,), ps);
            @test all(θs .<= vc.Θ_SAT);

            k1 = SH.relative_soil_k.((vc,), θs);
            k2 = SH.relative_soil_k.((vc,), true, ps);
            @test all(k1 .<= 1);
            @test all(k2 .<= 1);
            @test all(k1 .≈ k2);
        end;
    end;

    #=
    @testset "Vulnerability curve (fitting)" begin
        for vg in [NS.VanGenuchten{Float64}("Loam"), NS.VanGenuchten{Float64}("Sand")]
            bc = NS.BrooksCorey{Float64}(vg);
            @test true;
        end;
    end;
    =#

    @testset "Trace gas diffusion" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);

        # set the soil to be not saturated
        SPAC.update!(config, spac; swcs = (0.3, 0.3, 0.3, 0.3, 0.3));
        spac.SOILS[5].state.ns[4] = spac.AIRS[1].state.p_air * (spac.SOILS[5].state.vc.Θ_SAT - spac.SOILS[5].state.θ) * spac.SOILS[5].auxil.δz / (CS.GAS_R() * spac.SOILS[5].auxil.t);
        spac.SOILS[5].state.ns[5] = spac.AIRS[1].state.p_air * (spac.SOILS[5].state.vc.Θ_SAT - spac.SOILS[5].state.θ) * spac.SOILS[5].auxil.δz / (CS.GAS_R() * spac.SOILS[5].auxil.t);
        SH.volume_balance!(config, spac);
        SPAC.update_substep_auxils!(spac);
        SH.trace_gas_diffusion!(config, spac);

        # total mole into the layer should be the same as the derivative of the total mole of all layers
        for j in 1:5
            @test -spac.SOIL_BULK.auxil.dndt[1,j] ≈ sum([soil.auxil.∂n∂t[j] for soil in spac.SOILS]);
        end;

        # set one soil layer to be oversaturated
        spac.SOILS[3].state.θ = spac.SOILS[3].state.vc.Θ_SAT;
        SPAC.update_substep_auxils!(spac);
        SH.trace_gas_diffusion!(config, spac);

        # total mole into the layer should be the same as the derivative of the total mole of all layers
        for j in 1:5
            @test -spac.SOIL_BULK.auxil.dndt[1,j] ≈ sum([soil.auxil.∂n∂t[j] for soil in spac.SOILS]);
        end;

        # there is no dry air diffusion up to and from the third layer
        @test all(spac.SOIL_BULK.auxil.dndt[3,[1,2,4,5]] .== 0);
        @test all(spac.SOIL_BULK.auxil.dndt[4,[1,2,4,5]] .== 0);

        # set all soil layer to be oversaturated
        for j in 1:5
            spac.SOILS[j].state.θ = spac.SOILS[j].state.vc.Θ_SAT;
        end;
        SPAC.update_substep_auxils!(spac);
        SH.trace_gas_diffusion!(config, spac);

        # there is no dry air diffusion up to and from each layer
        for j in 1:5
            @test all(spac.SOIL_BULK.auxil.dndt[j,[1,2,4,5]] .== 0);
        end;
    end;

    @testset "Soil infiltration" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);

        # set the soil to be not saturated
        SPAC.update!(config, spac; swcs = (0.3, 0.3, 0.3, 0.3, 0.3));
        SPAC.update_substep_auxils!(spac);
        SH.soil_water_infiltration!(spac);

        # total water mole change by infiltration should be 0 (no water is lost or gained)
        @test sum([soil.auxil.∂θ∂t * soil.auxil.δz for soil in spac.SOILS]) ≈ 0 atol = 1e-10;

        # set the soil to be almost saturated (controlled by the soil water content)
        SPAC.update!(config, spac; swcs = (0.5, 0.5, 0.5, 0.5, 0.5));
        SPAC.update_substep_auxils!(spac);
        SH.soil_water_infiltration!(spac);

        # total water mole change by infiltration should be 0 (no water is lost or gained)
        @test sum([soil.auxil.∂θ∂t * soil.auxil.δz for soil in spac.SOILS]) ≈ 0 atol = 1e-10;
    end;

    @testset "Soil gas and water volume" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        SPAC.update!(config, spac; swcs = (0.3, 0.3, 0.3, 0.3, 0.3));
        SPAC.update_substep_auxils!(spac);

        # now the lowest layer to be short of dry air
        SH.trace_gas_diffusion!(config, spac);
        SH.volume_balance!(config, spac);
        for i in eachindex(spac.SOILS)
            soil = spac.SOILS[i];
            nmax = (spac.AIRS[1].state.p_air - PC.saturation_vapor_pressure(soil.auxil.t, soil.auxil.ψ * 1000000)) * soil.auxil.δz * (soil.state.vc.Θ_SAT - soil.state.θ) / (CS.GAS_R() * soil.auxil.t);
            ndry = soil.state.ns[1] + soil.state.ns[2] + soil.state.ns[4] + soil.state.ns[5];
            @test ndry <= nmax || ndry ≈ nmax;
        end;

        # set first layer to be saturated, and let the second layer to suck some water from the first layer
        spac.SOILS[1].state.θ = spac.SOILS[1].state.vc.Θ_SAT;
        spac.SOILS[1].state.ns .= 0;
        spac.SOILS[2].state.ns[5] = 0;
        θ_1 = spac.SOILS[1].state.θ;
        θ_2 = spac.SOILS[2].state.θ;
        SH.volume_balance!(config, spac);
        @test spac.SOILS[1].state.θ < θ_1;
        @test spac.SOILS[2].state.θ > θ_2;

        # set the last layer to be oversaturated in terms of dry air, and then air should move towards upper layers
        spac.SOILS[5].state.ns[4] = spac.AIRS[1].state.p_air * (spac.SOILS[5].state.vc.Θ_SAT - spac.SOILS[5].state.θ) * spac.SOILS[5].auxil.δz / (CS.GAS_R() * spac.SOILS[5].auxil.t);
        spac.SOILS[5].state.ns[5] = spac.AIRS[1].state.p_air * (spac.SOILS[5].state.vc.Θ_SAT - spac.SOILS[5].state.θ) * spac.SOILS[5].auxil.δz / (CS.GAS_R() * spac.SOILS[5].auxil.t);
        SH.volume_balance!(config, spac);
        for i in eachindex(spac.SOILS)
            soil = spac.SOILS[i];
            nmax = (spac.AIRS[1].state.p_air - PC.saturation_vapor_pressure(soil.auxil.t, soil.auxil.ψ * 1000000)) * soil.auxil.δz * (soil.state.vc.Θ_SAT - soil.state.θ) / (CS.GAS_R() * soil.auxil.t);
            ndry = soil.state.ns[1] + soil.state.ns[2] + soil.state.ns[4] + soil.state.ns[5];
            @test ndry <= nmax || ndry ≈ nmax;
        end;
    end;

    @testset "Surface runoff" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);

        # set every layer to be unsaturated
        SPAC.update!(config, spac; swcs = (0.3, 0.3, 0.3, 0.3, 0.3));
        SPAC.update_substep_auxils!(spac);
        SH.soil_water_runoff!(spac);
        @test spac.SOIL_BULK.auxil.runoff == 0;

        # set the second layer to be oversaturated (but not enough to cause runoff)
        SPAC.update!(config, spac; swcs = (0.3, spac.SOILS[2].state.vc.Θ_SAT + 0.01, 0.3, 0.3, 0.3));
        spac.SOILS[2].state.θ = spac.SOILS[2].state.vc.Θ_SAT + 0.01;
        SPAC.update_substep_auxils!(spac);
        SH.soil_water_runoff!(spac);
        @test spac.SOIL_BULK.auxil.runoff == 0;

        # set the second layer to be oversaturated (enough to cause runoff)
        SPAC.update!(config, spac; swcs = (0.3, 0.7, 0.3, 0.3, 0.3));
        spac.SOILS[2].state.θ = 0.7;
        SPAC.update_substep_auxils!(spac);
        SH.soil_water_runoff!(spac);
        runoffs = ((0.7 - spac.SOILS[2].state.vc.Θ_SAT) * spac.SOILS[2].auxil.δz + (0.3 - spac.SOILS[1].state.vc.Θ_SAT) * spac.SOILS[1].auxil.δz) * CS.ρ_H₂O(Float64) / CS.M_H₂O(Float64);
        @test spac.SOIL_BULK.auxil.runoff ≈ runoffs > 0;

        # set every layer to the oversaturated
        SPAC.update!(config, spac; swcs = (0.7, 0.7, 0.7, 0.7, 0.7));
        spac.SOILS[1].state.θ = 0.7;
        spac.SOILS[2].state.θ = 0.7;
        spac.SOILS[3].state.θ = 0.7;
        spac.SOILS[4].state.θ = 0.7;
        spac.SOILS[5].state.θ = 0.7;
        SPAC.update_substep_auxils!(spac);
        SH.soil_water_runoff!(spac);
        runoffs = sum([(0.7 - spac.SOILS[i].state.vc.Θ_SAT) * spac.SOILS[i].auxil.δz for i in 1:5]) * CS.ρ_H₂O(Float64) / CS.M_H₂O(Float64);
        @test spac.SOIL_BULK.auxil.runoff ≈ runoffs > 0;
    end;

    @testset "Water conservation" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        SPAC.update!(config, spac; swcs = (0.3, 0.3, 0.3, 0.3, 0.3));

        # the case of no root water uptake
        SPAC.update_substep_auxils!(spac);
        water_ini = sum([soil.state.θ * soil.auxil.δz for soil in spac.SOILS]) * CS.ρ_H₂O(Float64) / CS.M_H₂O(Float64) + sum([soil.state.ns[3] for soil in spac.SOILS]);
        for root in spac.ROOTS
            if typeof(root.xylem.auxil.flow) <: AbstractFloat
                root.xylem.auxil.flow = 0;
            else
                root.xylem.auxil.flow .= 0;
            end;
        end;
        SH.soil_profiles!(config, spac);
        SH.soil_budgets!(config, spac, 10.0);
        water_end = sum([soil.state.θ * soil.auxil.δz for soil in spac.SOILS]) * CS.ρ_H₂O(Float64) / CS.M_H₂O(Float64) + sum([soil.state.ns[3] for soil in spac.SOILS]);
        @test water_end + spac.SOIL_BULK.auxil.dndt[1,3] * 10 ≈ water_ini;

        # the case with root water uptake
        SPAC.update_substep_auxils!(spac);
        water_ini = sum([soil.state.θ * soil.auxil.δz for soil in spac.SOILS]) * CS.ρ_H₂O(Float64) / CS.M_H₂O(Float64) + sum([soil.state.ns[3] for soil in spac.SOILS]);
        for root in spac.ROOTS
            if typeof(root.xylem.auxil.flow) <: AbstractFloat
                root.xylem.auxil.flow = 1;
            else
                root.xylem.auxil.flow .= 1;
            end;
        end;
        SH.soil_profiles!(config, spac);
        SH.soil_budgets!(config, spac, 10.0);
        water_end = sum([soil.state.θ * soil.auxil.δz for soil in spac.SOILS]) * CS.ρ_H₂O(Float64) / CS.M_H₂O(Float64) + sum([soil.state.ns[3] for soil in spac.SOILS]);
        @test water_end + spac.SOIL_BULK.auxil.dndt[1,3] * 10 + length(spac.ROOTS) * 1 * 10 / spac.SOIL_BULK.state.area ≈ water_ini;

        # the case with root water uptake and surface runoff (water flow already defined above)
        spac.SOILS[1].state.θ = 0.7;
        SPAC.update_substep_auxils!(spac);
        water_ini = sum([soil.state.θ * soil.auxil.δz for soil in spac.SOILS]) * CS.ρ_H₂O(Float64) / CS.M_H₂O(Float64) + sum([soil.state.ns[3] for soil in spac.SOILS]);
        SH.soil_profiles!(config, spac);
        SH.soil_budgets!(config, spac, 10.0);
        water_end = sum([soil.state.θ * soil.auxil.δz for soil in spac.SOILS]) * CS.ρ_H₂O(Float64) / CS.M_H₂O(Float64) + sum([soil.state.ns[3] for soil in spac.SOILS]);
        @test water_end + spac.SOIL_BULK.auxil.dndt[1,3] * 10 + length(spac.ROOTS) * 1 * 10 / spac.SOIL_BULK.state.area + spac.SOIL_BULK.auxil.runoff ≈ water_ini;

        # the case with precipitation (water flow already defined above)
        SPAC.update_substep_auxils!(spac);
        spac.METEO.rain = 1;
        water_ini = sum([soil.state.θ * soil.auxil.δz for soil in spac.SOILS]) * CS.ρ_H₂O(Float64) / CS.M_H₂O(Float64) + sum([soil.state.ns[3] for soil in spac.SOILS]);
        SH.soil_profiles!(config, spac);
        SH.soil_budgets!(config, spac, 10.0);
        water_end = sum([soil.state.θ * soil.auxil.δz for soil in spac.SOILS]) * CS.ρ_H₂O(Float64) / CS.M_H₂O(Float64) + sum([soil.state.ns[3] for soil in spac.SOILS]);
        @test water_end + spac.SOIL_BULK.auxil.dndt[1,3] * 10 + length(spac.ROOTS) * 1 * 10 / spac.SOIL_BULK.state.area + spac.SOIL_BULK.auxil.runoff ≈ water_ini + spac.METEO.rain * 10;
    end;
end;
