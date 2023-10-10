@testset verbose = true "Tutorial" begin
    # Note that running these tutorials requires data to be prepared a priori, so we make an if statement to avoid it in automated CI tests.
    if isdir("/home/wyujie")
        include("tutorials/frontier.jl");
    end;

    include("tutorials/leaf_spectra.jl");
    include("tutorials/leaf_sif_spectra.jl");

    include("tutorials/spac_basic.jl");
    include("tutorials/spac_canopy_structure.jl");
    include("tutorials/spac_leaf_bio.jl");
    include("tutorials/spac_sun_sensor.jl");
end;
