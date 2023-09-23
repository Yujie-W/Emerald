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

    @testset "Pressure volume curves P(V)" begin
        for pv in [NS.ExponentialPVCurve{Float64}(), NS.LinearPVCurve{Float64}(), NS.SegmentedPVCurve{Float64}()]
            vs = collect(Float64, 0.21:0.01:1);
            ps = PH.capacitance_pressure.((pv,), vs, 298.15);

            @test all(-100 .<= ps);
        end;
    end;

    @testset "Pressure volume curves V(P)" begin
        for pv in [NS.ExponentialPVCurve{Float64}(), NS.LinearPVCurve{Float64}(), NS.SegmentedPVCurve{Float64}()]
            v1 = collect(Float64, 0.41:0.01:1);
            ps = PH.capacitance_pressure.((pv,), v1, 298.15);
            v2 = PH.capacitance_volume.((pv,), ps, 298.15);

            @test all(v2 .≈ v1);
        end;
    end;

    @testset "Read flow in/out of the xylem" begin
        config = NS.SPACConfiguration{Float64}();
        xylem = NS.XylemHydraulics(config);
        nssflow = NS.XylemHydraulicsAuxilNSS(config);
        ssflow = NS.XylemHydraulicsAuxilSS(config);

        @test PH.flow_in(xylem) == PH.flow_in(nssflow) == PH.flow_in(ssflow) == 0;
        @test PH.flow_out(xylem) == PH.flow_out(ssflow) == PH.flow_out(nssflow) == 0;

        PH.set_flow_profile!(xylem, 1.0);
        PH.set_flow_profile!(nssflow, 1.0);
        PH.set_flow_profile!(ssflow, 1.0);

        @test PH.flow_in(xylem) == PH.flow_in(nssflow) == PH.flow_in(ssflow) == 1;
        @test PH.flow_out(xylem) == PH.flow_out(ssflow) == PH.flow_out(nssflow) == 1;
    end;

end;
