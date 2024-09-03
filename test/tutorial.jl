using Test


@testset verbose = true "Tutorial" begin
    # Note that running these tutorials requires data to be prepared a priori, so we make an if statement to avoid it in automated CI tests.
    if isdir("/home/wyujie")
        include("tutorials/frontier.jl");
    end;

    # change the configurations
    include("tutorials/config/turn_off_refl_sif.jl");

    include("tutorials/leaf_quantum_yields.jl");
    include("tutorials/leaf_spectra.jl");
    include("tutorials/leaf_sif_spectra.jl");

    include("tutorials/plant_death.jl");

    include("tutorials/spac_basic.jl");
    include("tutorials/spac_canopy_structure.jl");
    include("tutorials/spac_leaf_bio.jl");
    include("tutorials/spac_par_wl.jl");
    include("tutorials/spac_sif.jl");
    include("tutorials/spac_stomatal_model.jl");
    include("tutorials/spac_sun_sensor.jl");
    include("tutorials/spac_sync_states.jl");
end;
