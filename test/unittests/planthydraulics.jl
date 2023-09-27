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

    @testset "Update pressure profile in xylem" begin
        config = NS.SPACConfiguration{Float64}();
        xylem = NS.XylemHydraulics(config);
        nssflow = NS.XylemHydraulicsAuxilNSS(config);
        ssflow = NS.XylemHydraulicsAuxilSS(config);
        PH.set_flow_profile!(xylem, 1.0);
        PH.set_flow_profile!(nssflow, 1.0);
        PH.set_flow_profile!(ssflow, 1.0);

        PH.xylem_pressure_profile!(xylem, 298.15);
        PH.xylem_pressure_profile!(xylem.state, nssflow, 298.15);
        PH.xylem_pressure_profile!(xylem.state, ssflow, 298.15);

        @test xylem.auxil.pressure[end] == nssflow.pressure[end] == ssflow.pressure[end] < 0;
        @test all(xylem.auxil.pressure[2:end] .< 0);
        @test all(nssflow.pressure[2:end] .< 0);
        @test all(ssflow.pressure[2:end] .< 0);
    end;

    @testset "Root flow and pressure profiles" begin
        config = NS.SPACConfiguration{Float64}();
        root = NS.Root(config);
        soil = NS.SoilLayer{Float64}();
        flow = 1.0;
        PH.set_flow_profile!(root.xylem, flow);
        PH.root_pressure_profile!(root, soil);

        p_target = root.xylem.auxil.pressure[end];
        PH.root_flow_profile!(root, soil, p_target);
        f_target = PH.flow_out(root.xylem);
        PH.root_pressure_profile!(root, soil);

        @test PH.flow_out(root.xylem) ≈ f_target;
        @test root.xylem.auxil.pressure[end] ≈ p_target;
    end;

    @testset "Stem flow and pressure profiles" begin
        config = NS.SPACConfiguration{Float64}();
        stem = NS.Stem(config);
        flow = 1.0;
        PH.set_flow_profile!(stem.xylem, flow);
        PH.stem_pressure_profile!(stem, -0.1);

        @test PH.flow_out(stem.xylem) == flow;
        @test all(stem.xylem.auxil.pressure .<= -0.1);
    end;

    @testset "Leaf flow and pressure profile" begin
        config1 = NS.SPACConfiguration{Float64}();
        config2 = NS.SPACConfiguration{Float64}(STEADY_STATE_FLOW = false);
        leaf1 = NS.Leaf(config);
        leaf2 = NS.Leaf(config);
        PH.leaf_pressure_profile!(leaf1, -0.1);
        PH.leaf_pressure_profile!(leaf2, -0.1);

        @test all(leaf1.xylem.auxil.pressure .<= -0.1);
        @test all(leaf2.xylem.auxil.pressure .<= -0.1);
    end;

end;
