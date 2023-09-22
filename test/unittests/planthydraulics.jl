using Test
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.PlantHydraulics as PH


@testset verbose = true "Plant Hydraulics Model" begin
    @testset "Vulnerability curves K(P)" begin
        for vc in [NS.ComplexVC{Float64}(), NS.LogisticVC{Float64}(), NS.PowerVC{Float64}(), NS.WeibullVC{Float64}()]
            ps = collect(Float64, 0:-0.1:-10);
            ks = PH.relative_xylem_k.((vc,), ps);

            @test all(0 .<= ks .<= 1);
            @test ks[1] == 1;
        end;
    end;

    @testset "Vulnerability curves P(K)" begin
        for vc in [NS.ComplexVC{Float64}(), NS.LogisticVC{Float64}(), NS.PowerVC{Float64}(), NS.WeibullVC{Float64}()]
            p1 = collect(Float64, 0:-0.1:-5);
            ks = PH.relative_xylem_k.((vc,), p1);
            p2 = PH.xylem_pressure.((vc,), ks);

            @test all(p2 .≈ p1);
        end;
    end;

end;
