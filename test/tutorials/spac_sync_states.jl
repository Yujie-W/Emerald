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
    #     1. LAI > 0 and SAI > 0 (vegetation with leaves, e.g. forest)
    #     2. LAI = 0 and SAI > 0 (vegetation without leaves, e.g. winter or dry season)
    #     2. LAI > 0 and SAI = 0 (vegetation with leaves only, e.g. grassland or crops)
    #     3. LAI = 0 and SAI = 0 (bare soil)
    # In the following example, we will run the steps to see if our SPAC model can handle these scenarios.
    # create a SPAC to work on

    @testset "Run SPAC with sync_state!" begin
        config = EmeraldLand.Namespace.SPACConfiguration(FT);
        spac = EmeraldLand.Namespace.BulkSPAC(config);
        spac_bak = deepcopy(spac);
        EmeraldLand.SPAC.initialize_spac!(config, spac);
        state_0 = EmeraldLand.Namespace.BulkSPACStates(spac);

        # compare the models using different modes (regular or state sync modes)
        state_1 = deepcopy(state_0);
        state_2 = deepcopy(state_0);
        spac_1 = deepcopy(spac);
        for i in 1:3
            spac_2 = deepcopy(spac_bak);
            EmeraldLand.SPAC.initialize_spac!(config, spac_2, state_2);
            EmeraldLand.SPAC.spac!(config, spac_1, FT(3600));
            EmeraldLand.SPAC.spac!(config, spac_2, FT(3600));
            EmeraldLand.Namespace.sync_state!(spac_1, state_1);
            EmeraldLand.Namespace.sync_state!(spac_2, state_2);
            EmeraldUtility.StructEqual.compare_struct!(spac_1, spac_2; approximation = false, show_diff_msg = true);
            @test EmeraldUtility.StructEqual.compare_struct!(spac_1, spac_2; show_diff_msg = false) == 0;
        end;
    end;

    @testset "SPAC with LAI > 0 and SAI = 0" begin
        config = EmeraldLand.Namespace.SPACConfiguration(FT);
        spac = EmeraldLand.Namespace.BulkSPAC(config);
        EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 2.0, sai = 0.0);
        spac_bak = deepcopy(spac);
        EmeraldLand.SPAC.initialize_spac!(config, spac);
        state_0 = EmeraldLand.Namespace.BulkSPACStates(spac);

        # compare the models using different modes (regular or state sync modes)
        state_1 = deepcopy(state_0);
        state_2 = deepcopy(state_0);
        spac_1 = deepcopy(spac);
        spac_2 = deepcopy(spac_bak);
        EmeraldLand.SPAC.initialize_spac!(config, spac_2, state_2);
        EmeraldLand.SPAC.spac!(config, spac_1, FT(3600));
        EmeraldLand.SPAC.spac!(config, spac_2, FT(3600));
        EmeraldLand.Namespace.sync_state!(spac_1, state_1);
        EmeraldLand.Namespace.sync_state!(spac_2, state_2);
        @test EmeraldUtility.StructEqual.compare_struct!(spac_1, spac_2; show_diff_msg = false) == 0;
    end;

    @testset "SPAC with LAI = 0 and SAI > 0" begin
        config = EmeraldLand.Namespace.SPACConfiguration(FT);
        spac = EmeraldLand.Namespace.BulkSPAC(config);
        EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 0.0, sai = 0.5);
        spac_bak = deepcopy(spac);
        EmeraldLand.SPAC.initialize_spac!(config, spac);
        state_0 = EmeraldLand.Namespace.BulkSPACStates(spac);

        # compare the models using different modes (regular or state sync modes)
        state_1 = deepcopy(state_0);
        state_2 = deepcopy(state_0);
        spac_1 = deepcopy(spac);
        spac_2 = deepcopy(spac_bak);
        EmeraldLand.SPAC.initialize_spac!(config, spac_2, state_2);
        EmeraldLand.SPAC.spac!(config, spac_1, FT(3600));
        EmeraldLand.SPAC.spac!(config, spac_2, FT(3600));
        EmeraldLand.Namespace.sync_state!(spac_1, state_1);
        EmeraldLand.Namespace.sync_state!(spac_2, state_2);
        @test EmeraldUtility.StructEqual.compare_struct!(spac_1, spac_2; show_diff_msg = false) == 0;
    end;

    @testset "SPAC with LAI = 0 and SAI = 0" begin
        config = EmeraldLand.Namespace.SPACConfiguration(FT);
        spac = EmeraldLand.Namespace.BulkSPAC(config);
        EmeraldLand.SPAC.prescribe_traits!(config, spac; lai = 0.0, sai = 0.0);
        spac_bak = deepcopy(spac);
        EmeraldLand.SPAC.initialize_spac!(config, spac);
        state_0 = EmeraldLand.Namespace.BulkSPACStates(spac);

        # compare the models using different modes (regular or state sync modes)
        state_1 = deepcopy(state_0);
        state_2 = deepcopy(state_0);
        spac_1 = deepcopy(spac);
        spac_2 = deepcopy(spac_bak);
        EmeraldLand.SPAC.initialize_spac!(config, spac_2, state_2);
        EmeraldLand.SPAC.spac!(config, spac_1, FT(3600));
        EmeraldLand.SPAC.spac!(config, spac_2, FT(3600));
        EmeraldLand.Namespace.sync_state!(spac_1, state_1);
        EmeraldLand.Namespace.sync_state!(spac_2, state_2);
        @test EmeraldUtility.StructEqual.compare_struct!(spac_1, spac_2; show_diff_msg = false) == 0;
    end;
end;
