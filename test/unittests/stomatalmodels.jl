using Test
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.Photosynthesis as PS
import Emerald.EmeraldLand.StomatalModels as SM


@testset verbose = true "StomatalModels.jl" begin
    @testset "Empirical equations" begin
        config = NS.SPACConfiguration{Float64}();
        leaf = NS.Leaf(config);
        air = NS.AirLayer{Float64}();
        leaf.flux.auxil.ppar_sunlit .= 100.0;
        leaf.flux.auxil.ppar_shaded = 100.0;
        leaf.flux.auxil.g_CO₂_shaded = 0.02;
        leaf.flux.auxil.g_CO₂_sunlit .= 0.02;
        PS.leaf_photosynthesis!(leaf, air, NS.GCO₂Mode(), 1.0; rd_only = false);

        for sm in [NS.BallBerrySM{Float64}(), NS.GentineSM{Float64}(), NS.LeuningSM{Float64}(), NS.MedlynSM{Float64}()]
            leaf.flux.state.stomatal_model = sm;
            gsh = SM.empirical_equation(sm, leaf, air);
            gsl = SM.empirical_equation(sm, leaf, air, 1);
            @test gsh > 0.0;
            @test gsl > 0.0;
        end;
    end;

end;
