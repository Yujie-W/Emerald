using Emerald;
using Test;


@testset "Run SPAC Response to SWC" begin
    FT = Float64;

    # To run SPAC model, users need to define the following to proceed, and then initialize and run the model:
    #     - Configuration struct
    #     - SPAC struct
    # For details about how to modify the SPAC struct or configurations, please refer to other tutorials.
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, FT(1));
    @test true;

    # Change SWC and the junction water content accordingly
    θ_min = spac.soils[1].trait.vc.Θ_RES;
    θ_max = spac.soils[1].trait.vc.Θ_SAT;
    for θ in 0.09:-0.0001:0.085
        for s in spac.soils
            s.state.θ = θ;
        end;
        psoil = EmeraldLand.SoilHydraulics.soil_ψ_25(spac.soils[1].trait.vc, θ) * EmeraldLand.PhysicalChemistry.relative_surface_tension(spac.soils[1].s_aux.t);
        spac.plant.junction.state.v_storage = EmeraldLand.PlantHydraulics.capacitance_volume(spac.plant.junction.trait.pv, psoil, spac.plant.junction.s_aux.t) * spac.plant.junction.trait.v_max;
        EmeraldLand.SPAC.initialize_spac!(config, spac);
        EmeraldLand.SPAC.spac!(config, spac, 3600);
        @show θ EmeraldLand.SPAC.GPP(spac);
    end;
end;
