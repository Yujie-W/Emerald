using Test
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.PlantHydraulics as PH
import Emerald.EmeraldLand.SPAC


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
        junc = NS.JunctionCapacitor{Float64}();
        flow = 1.0;
        PH.set_flow_profile!(root.xylem, flow);
        PH.root_pressure_profile!(soil, root, junc);

        p_target = root.xylem.auxil.pressure[end];
        junc.auxil.pressure = p_target;
        PH.root_flow_profile!(config, root, soil, junc);
        f_target = PH.flow_out(root);
        PH.root_pressure_profile!(soil, root, junc);

        @test PH.flow_out(root) ≈ f_target;
        @test root.xylem.auxil.pressure[end] ≈ p_target;
    end;

    @testset "Stem flow and pressure profiles" begin
        config = NS.SPACConfiguration{Float64}();
        stem = NS.Stem(config);
        flow = 1.0;
        PH.set_flow_profile!(stem.xylem, flow);
        PH.stem_pressure_profile!(stem, -0.1);

        @test PH.flow_out(stem) == flow;
        @test all(stem.xylem.auxil.pressure .<= -0.1);
    end;

    @testset "Leaf flow and pressure profile" begin
        config1 = NS.SPACConfiguration{Float64}();
        config2 = NS.SPACConfiguration{Float64}(STEADY_STATE_FLOW = false);
        leaf1 = NS.Leaf(config1);
        leaf2 = NS.Leaf(config2);
        PH.leaf_pressure_profile!(config1, leaf1, -0.1);
        PH.leaf_pressure_profile!(config2, leaf2, -0.1);

        @test all(leaf1.xylem.auxil.pressure .<= -0.1);
        @test all(leaf2.xylem.auxil.pressure .<= -0.1);
    end;

    @testset "Plant hydraulics (steady state)" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        for leaf in spac.LEAVES
            leaf.flux.state.g_H₂O_s_sunlit .= 0.1;
            leaf.flux.state.g_H₂O_s_shaded = 0.1;
        end;
        PH.plant_water_budget!(spac, 1.0);
        PH.plant_flow_profile!(config, spac);
        PH.plant_pressure_profile!(config, spac);

        # make sure the flow in and out are the same
        for i in eachindex(spac.LEAVES)
            @test PH.flow_out(spac.LEAVES[i]) == PH.flow_in(spac.LEAVES[i]) == PH.flow_out(spac.BRANCHES[i]) == PH.flow_in(spac.BRANCHES[i]);
        end;
        @test PH.flow_in(spac.TRUNK) == PH.flow_out(spac.TRUNK) == sum([PH.flow_in(branch) for branch in spac.BRANCHES]);
        for i in eachindex(spac.ROOTS)
            @test PH.flow_in(spac.ROOTS[i]) == PH.flow_out(spac.ROOTS[i]);
        end;

        # make sure the water in the junction capacitor changes with the total water in (from roots) and total water out (to air)
        Σf_root = sum([PH.flow_in(root) for root in spac.ROOTS]);
        Σf_leaf = sum([PH.flow_out(leaf) + leaf.capacitor.auxil.flow for leaf in spac.LEAVES]);
        q1_junc = spac.JUNCTION.state.v_storage;
        PH.plant_water_budget!(spac, 1.0);
        q2_junc = spac.JUNCTION.state.v_storage;
        @test q2_junc - q1_junc ≈ Σf_root - Σf_leaf;
    end;

    @testset "Plant hydraulics (non-steady state)" begin
        config = NS.SPACConfiguration{Float64}(STEADY_STATE_FLOW = false);
        spac = NS.MultiLayerSPAC(config);
        SPAC.initialize!(config, spac);
        for leaf in spac.LEAVES
            leaf.flux.state.g_H₂O_s_sunlit .= 0.1;
            leaf.flux.state.g_H₂O_s_shaded = 0.1;
        end;
        PH.plant_water_budget!(spac, 1.0);
        PH.plant_flow_profile!(config, spac);
        PH.plant_pressure_profile!(config, spac);

        # make sure the flow in and out are the same
        for i in eachindex(spac.LEAVES)
            @test PH.flow_out(spac.BRANCHES[i]) == PH.flow_in(spac.LEAVES[i]);
        end;
        @test PH.flow_out(spac.TRUNK) == sum([PH.flow_in(branch) for branch in spac.BRANCHES]);

        # make sure the water in the leaf capacitor changes with the total water in (from stem) and total water out (to air)
        q1s = [leaf.capacitor.state.v_storage for leaf in spac.LEAVES];
        fis = [PH.flow_in(leaf) for leaf in spac.LEAVES];
        fos = [PH.flow_out(leaf) for leaf in spac.LEAVES];
        PH.plant_water_budget!(spac, 1.0);
        PH.plant_flow_profile!(config, spac);
        PH.plant_pressure_profile!(config, spac);
        q2s = [leaf.capacitor.state.v_storage for leaf in spac.LEAVES];
        @test all(q2s .- q1s .≈ fis .- fos);

        # make sure the water in the branch capacitor changes with the total water in (from trunk) and total water out (to leaf)
        q1s = [sum(branch.xylem.state.v_storage) for branch in spac.BRANCHES];
        fis = [PH.flow_in(branch) for branch in spac.BRANCHES];
        fos = [PH.flow_out(branch) for branch in spac.BRANCHES];
        PH.plant_water_budget!(spac, 1.0);
        PH.plant_flow_profile!(config, spac);
        PH.plant_pressure_profile!(config, spac);
        q2s = [sum(branch.xylem.state.v_storage) for branch in spac.BRANCHES];
        @test all(q2s .- q1s .≈ fis .- fos);

        # make sure the water in the trunk capacitor changes with the total water in (from junction) and total water out (to branches)
        q1_trunk = sum(spac.TRUNK.xylem.state.v_storage);
        f_junc = PH.flow_in(spac.TRUNK);
        f_stem = PH.flow_out(spac.TRUNK);
        PH.plant_water_budget!(spac, 1.0);
        PH.plant_flow_profile!(config, spac);
        PH.plant_pressure_profile!(config, spac);
        q2_trunk = sum(spac.TRUNK.xylem.state.v_storage);
        @test q2_trunk - q1_trunk ≈ f_junc - f_stem;

        # make sure the water in the junction capacitor changes with the total water in (from roots) and total water out (to trunk)
        q1_junc = spac.JUNCTION.state.v_storage;
        Σf_root = sum([PH.flow_out(root) for root in spac.ROOTS]);
        Σf_stem = PH.flow_in(spac.TRUNK);
        PH.plant_water_budget!(spac, 1.0);
        PH.plant_flow_profile!(config, spac);
        PH.plant_pressure_profile!(config, spac);
        q2_junc = spac.JUNCTION.state.v_storage;
        @test q2_junc - q1_junc ≈ Σf_root - Σf_stem;

        # make sure the water in the root capacitor changes with the total water in (from soil) and total water out (to junction)
        q1s = [sum(root.xylem.state.v_storage) for root in spac.ROOTS];
        fis = [PH.flow_in(root) for root in spac.ROOTS];
        fos = [PH.flow_out(root) for root in spac.ROOTS];
        PH.plant_water_budget!(spac, 1.0);
        PH.plant_flow_profile!(config, spac);
        PH.plant_pressure_profile!(config, spac);
        q2s = [sum(root.xylem.state.v_storage) for root in spac.ROOTS];
        @test all(q2s .- q1s .≈ fis .- fos);
    end;
end;
