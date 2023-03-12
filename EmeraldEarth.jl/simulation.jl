#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to run the SPAC
#
#######################################################################################################################################################################################################
"""

    simulation!(mat_spac::Matrix; threads::Int = 12)
    simulation!(spac::Union{Nothing,MonoMLTreeSPAC})

Run simulations on SPAC, given
- `mat_spac` Matrix of SPAC
- `threads` Number of threadings
- `spac` SPAC or nothing

---
# Example
```julia
using Emerald;
dts = EmeraldEarth.LandDatasets{Float64}("gm2", 2020);
mat = EmeraldEarth.spac_grids(dts);
wd = EmeraldEarth.ERA5SingleLevelsDriver();
@time EmeraldEarth.prescribe!(mat, dts, wd, 1);
@time EmeraldEarth.simulation!(mat);
```

"""
function simulation! end

simulation!(mat_spac::Matrix; threads::Int = 12) = (
    add_threads!(threads);

    @tinfo "Running the global simulations in multiple threads...";
    @showprogress pmap(simulation!, mat_spac);

    return nothing
);

simulation!(spac::Union{Nothing,MonoMLTreeSPAC}) = (
    for _i in 1:30
        soil_plant_air_continuum!(spac, 120; p_on = false, t_on = false, θ_on = false);
    end;

    return nothing
);
