#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to run the SPAC
#
#######################################################################################################################################################################################################
"""

    simulation!(mat_spac::Matrix)
    simulation!(spac::Union{Nothing,MonoMLTreeSPAC})

Run simulations on SPAC, given
- `mat_spac` Matrix of SPAC
- `spac` SPAC or nothing

---
# Example
```julia
using Emerald;
dts = EmeraldEarth.LandDatasets{Float64}("gm2",2020);
mat = EmeraldEarth.spac_grids(dts);
wd = EmeraldEarth.ERA5SingleLevelsDriver();
@time EmeraldEarth.prescribe!(mat,dts,wd,1);
@time EmeraldEarth.global_simulation!(mat);
```

"""
function simulation! end

simulation!(mat_spac::Matrix) = (
    for _i in eachindex(mat_spac)
        _spac = mat_spac[_i];
        if _spac isa MonoMLTreeSPAC
            simulation!(_spac);
        end;
    end;
    #run_time_step!.(mat_spac);

    return nothing
);

simulation!(spac::Union{Nothing,MonoMLTreeSPAC}) = (
    for _i in 1:30
        soil_plant_air_continuum!(spac, 120; p_on = false, t_on = false, θ_on = false);
    end;

    return nothing
);
