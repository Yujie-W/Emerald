using Test
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.SPAC as SPAC


@testset verbose = true "Namespace.jl" begin
    @testset "State constructors" begin
        config = NS.SPACConfiguration{Float64}();
        spac = NS.BulkSPAC(config);
        SPAC.initialize_states!(config, spac);
        SPAC.initialize_spac!(config, spac);
        NS.t_aux!(config, spac.canopy);
        NS.s_aux!(config, spac.canopy);
        states = NS.BulkSPACStates(spac);
        @test true;
        NS.sync_state!(spac, states);
        @test true;
        NS.sync_state!(states, spac);
        @test true;
    end;
end;
