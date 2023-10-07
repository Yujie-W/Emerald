using Test
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.SoilHydraulics as SH
import Emerald.EmeraldLand.SPAC

@testset verbose = true "Soil Hydraulics Model" begin
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

    @testset "Vulnerability curve (fitting)" begin
        for vg in [NS.VanGenuchten{Float64}("Loam"), NS.VanGenuchten{Float64}("Sand")]
            bc = NS.BrooksCorey{Float64}(vg);
            @test true;
        end;
    end;

    @testset "Trace Gas Diffusion" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        SPAC.update!(config, spac; swcs = (0.3, 0.3, 0.3, 0.3, 0.3));
        SPAC.update_substep_auxils!(spac);
        SH.trace_gas_diffusion!(config, spac);

        # there is not NaN in the derivative of mole
        @test all(!isnan, spac.SOIL_BULK.auxil.dndt);
        for soil in spac.SOILS
            @test all(!isnan, soil.auxil.∂n∂t);
        end;

        # total mole into the layer should be the same as the derivative of the total mole of all layers
        for j in 1:5
            @test -spac.SOIL_BULK.auxil.dndt[1,j] ≈ sum([soil.auxil.∂n∂t[j] for soil in spac.SOILS]);
        end;
    end;

end;
