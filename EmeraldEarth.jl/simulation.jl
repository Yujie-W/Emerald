#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to run the SPAC
#
#######################################################################################################################################################################################################
"""

    simulation!(mat_spac::Matrix{Union{Nothing,MonoMLTreeSPAC{FT}}}) where {FT<:AbstractFloat}
    simulation!(spac::Union{Nothing,MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat}

Run simulations on SPAC, given
- `mat_spac` Matrix of SPAC
- `spac` SPAC or nothing

---
# Example
```julia
using Emerald;

FT = Float64;
@time EmeraldEarth.add_threads!(4, FT);
@time dts = EmeraldEarth.LandDatasets{FT}("gm2", 2020);
@time mat = EmeraldEarth.spac_grids(dts);
@time wdr = EmeraldEarth.ERA5SingleLevelsDriver();
@time wdx = EmeraldEarth.weather_grids(dts, wdr, 6);

@time EmeraldEarth.prescribe!(mat, dts, wdr, 6);
@time mat = EmeraldEarth.simulation!(mat; threads = 120);
```

"""
function simulation! end

simulation!(mat_spac::Matrix{Union{Nothing,MonoMLTreeSPAC{FT}}}) where {FT<:AbstractFloat} = (
    @tinfo "Running the global simulations in multiple threads...";
    _spacs = @showprogress pmap(simulation!, mat_spac);

    return _spacs
);

simulation!(spac::Union{Nothing,MonoMLTreeSPAC{FT}}) where {FT<:AbstractFloat} = (
    for _i in 1:10
        soil_plant_air_continuum!(spac, 360; p_on = false, t_on = false, θ_on = false);
    end;

    return spac
);
