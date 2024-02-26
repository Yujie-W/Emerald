using Emerald;
using Test;


@testset verbose = true "Sync State and Run SPAC" begin
    FT = Float64;

    # The SPAC sync state functon is meant to sync the states of the SPAC system. This feature is meant for global simulations as the cached SPAC struct could be different among grids.
    # The regular steps for any grid are
    #     1. Initialize the SPAC struct and states
    #     2. Run the SPAC model and return the states
    #     3. Sync the states back to the SPAC struct
    #     4. Run the SPAC model again using the new states in the next time step
    # Currently, there are a few scenarios based on LAI and SAI:
    #     1. LAI > 0 and SAI > 0 (vegetation with leaves)
    #     2. LAI = 0 and SAI > 0 (vegetation without leaves, e.g. winter or dry season)
    #     3. LAI = 0 and SAI = 0 (bare soil)
    # In the following example, we will run the steps to see if our SPAC model can handle these scenarios.
    # create a SPAC to work on

    @testset "SPAC with LAI > 0 and SAI > 0" begin
        config = EmeraldLand.Namespace.SPACConfiguration{FT}();
        spac = EmeraldLand.Namespace.BulkSPAC(config);
        EmeraldLand.SPAC.initialize_states!(config, spac);
        EmeraldLand.SPAC.initialize_spac!(config, spac);
        EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 3.0, sai = 0.5);
        @test true;

        # create a state to sync to
        state = EmeraldLand.Namespace.BulkSPACStates(spac);
        @test true;

        # run the spac model for 1 hour
        spac_1 = deepcopy(spac);
        EmeraldLand.SPAC.spac!(config, spac_1, FT(3600));
        @test true;

        # sync the states back to the original spac and run the model for 1 hour
        spac_2 = deepcopy(spac);
        EmeraldLand.Namespace.sync_state!(state, spac_2);
        EmeraldLand.SPAC.spac!(config, spac_2, FT(3600));
        @test true;
    end;

    @testset "SPAC with LAI > 0 and SAI > 0" begin
        config = EmeraldLand.Namespace.SPACConfiguration{FT}();
        spac = EmeraldLand.Namespace.BulkSPAC(config);
        EmeraldLand.SPAC.initialize_states!(config, spac);
        EmeraldLand.SPAC.initialize_spac!(config, spac);
        EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 0.0, sai = 0.5);
        @test true;

        # create a state to sync to
        state = EmeraldLand.Namespace.BulkSPACStates(spac);
        @test true;

        # run the spac model for 1 hour
        spac_1 = deepcopy(spac);
        EmeraldLand.SPAC.spac!(config, spac_1, FT(3600));
        @test true;

        # sync the states back to the original spac and run the model for 1 hour
        spac_2 = deepcopy(spac);
        EmeraldLand.Namespace.sync_state!(state, spac_2);
        EmeraldLand.SPAC.spac!(config, spac_2, FT(3600));
        @test true;
    end;

    @testset "SPAC with LAI > 0 and SAI = 0" begin
        config = EmeraldLand.Namespace.SPACConfiguration{FT}();
        spac = EmeraldLand.Namespace.BulkSPAC(config);
        EmeraldLand.SPAC.initialize_states!(config, spac);
        EmeraldLand.SPAC.initialize_spac!(config, spac);
        EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 0.0, sai = 0.0);
        @test true;

        # create a state to sync to
        state = EmeraldLand.Namespace.BulkSPACStates(spac);
        @test true;

        # run the spac model for 1 hour
        spac_1 = deepcopy(spac);
        EmeraldLand.SPAC.spac!(config, spac_1, FT(3600));
        @test true;

        # sync the states back to the original spac and run the model for 1 hour
        spac_2 = deepcopy(spac);
        EmeraldLand.Namespace.sync_state!(state, spac_2);
        EmeraldLand.SPAC.spac!(config, spac_2, FT(3600));
        @test true;
    end;

    @testset "LAI and SAI changed (> 0)" begin
        config = EmeraldLand.Namespace.SPACConfiguration{FT}();
        spac = EmeraldLand.Namespace.BulkSPAC(config);
        EmeraldLand.SPAC.initialize_states!(config, spac);
        EmeraldLand.SPAC.initialize_spac!(config, spac);
        EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 0.0, sai = 0.0);
        @test true;

        # create a state to sync to
        state = EmeraldLand.Namespace.BulkSPACStates(spac);
        @test true;

        # run the spac model for 1 hour
        spac_1 = deepcopy(spac);
        EmeraldLand.SPAC.prescribe_traits!(config, spac_1; lai = 3.0, sai = 0.5);
        EmeraldLand.SPAC.spac!(config, spac_1, FT(3600));
        @test true;

        # sync the states back to the original spac and run the model for 1 hour
        spac_2 = deepcopy(spac);
        EmeraldLand.Namespace.sync_state!(state, spac_2);
        EmeraldLand.SPAC.prescribe_traits!(config, spac_2; lai = 3.0, sai = 0.5);
        EmeraldLand.SPAC.spac!(config, spac_2, FT(3600));
        @test true;
    end;

    @testset "LAI and SAI changed (= 0)" begin
        config = EmeraldLand.Namespace.SPACConfiguration{FT}();
        spac = EmeraldLand.Namespace.BulkSPAC(config);
        EmeraldLand.SPAC.initialize_states!(config, spac);
        EmeraldLand.SPAC.initialize_spac!(config, spac);
        EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 3.0, sai = 0.5);
        @test true;

        # create a state to sync to
        state = EmeraldLand.Namespace.BulkSPACStates(spac);
        @test true;

        # run the spac model for 1 hour
        spac_1 = deepcopy(spac);
        EmeraldLand.SPAC.prescribe_traits!(config, spac_1; lai = 0.0, sai = 0.0);
        EmeraldLand.SPAC.spac!(config, spac_1, FT(3600));
        @test true;

        # sync the states back to the original spac and run the model for 1 hour
        spac_2 = deepcopy(spac);
        EmeraldLand.Namespace.sync_state!(state, spac_2);
        EmeraldLand.SPAC.prescribe_traits!(config, spac_2; lai = 0.0, sai = 0.0);
        EmeraldLand.SPAC.spac!(config, spac_2, FT(3600));
        @test true;
    end;

    @test true;
end;
