using Emerald
using Test

@testset "Change Xylem Element Number" begin
    FT = Float64;

    for N in [1; collect(5:5:200)]
        config = EmeraldLand.Namespace.SPACConfiguration(FT);
        config.DIM_XYLEM = N;
        stem = EmeraldLand.Namespace.Stem(config);
        stem.xylem.auxil.flow = 2.2;
        stem.xylem.auxil.e_crit = EmeraldLand.PlantHydraulics.critical_flow(config, stem.xylem, stem.energy.s_aux.t);
        EmeraldLand.PlantHydraulics.stem_pressure_profile!(stem, FT(-2));
        @test !isnan(stem.xylem.auxil.pressure[end]);
    end;
end;
