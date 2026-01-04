using PkgUtility.DataIO: read_csv
using Test

import Emerald.Namespace as ENS
import Emerald.ResearchTools as RTS


@testset "Emerald ResearchTools" verbose = true begin
    df3 = read_csv(joinpath(@__DIR__, "../..", "data/examples", "C3-ACi.csv"));
    df4 = read_csv(joinpath(@__DIR__, "../..", "data/examples", "C4-ACi.csv"));
    df3.T_LEAF .+= 273.15;  # convert to Kelvin
    df4.T_LEAF .+= 273.15;  # convert to Kelvin

    @testset "C3" begin
        config = ENS.SPACConfig(Float64);
        result = RTS.ACi.aci_fit!(config, df3, "C3", ["Vcmax25", "Jmax25", "Γstar25", "Rd25"]);
    end;
end;
