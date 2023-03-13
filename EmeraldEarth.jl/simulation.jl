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
@time EmeraldEarth.add_threads!(20, FT);
@time dts = EmeraldEarth.LandDatasets{FT}("gm2", 2020);
@time mat = EmeraldEarth.gm_grids(dts);
@time wdr = EmeraldEarth.ERA5SingleLevelsDriver();
@time wdx = EmeraldEarth.wd_grids(dts, wdr, 6);
@time sts = EmeraldEarth.simulation!(mat, wdx);

nansts = zeros(Bool, size(sts));
for i in eachindex(sts)
    if !isnothing(sts[i]) && isnan(sts[i])
        nansts[i] = true;
    end;
end;
@show sum(nansts);
sts[isnothing.(sts)] .= NaN;
heatmap(sts')

EmeraldEarth.simulation!(mat[339,36], wdx[339,36])
EmeraldEarth.simulation!(mat[296,120], wdx[296,120])
```

"""
function simulation! end

simulation!(gm_mat::Matrix{Union{Nothing,Dict{String,Any}}}, wd_mat::Matrix{Union{Nothing,Dict{String,Any}}}) = (
    @tinfo "Running the global simulations in multiple threads...";
    _states = @showprogress pmap(simulation!, gm_mat, wd_mat);

    return _states
);

simulation!(gm_params::Nothing, wd_params::Nothing) = nothing;

simulation!(gm_params::Dict{String,Any}, wd_params::Dict{String,Any}) = (
    synchronize_cache!(gm_params, wd_params);
    for _i in 1:10
        soil_plant_air_continuum!(CACHE_SPAC, 360; p_on = false, t_on = false, θ_on = false);
    end;

    # if wd_params["RAD_DIR"] > 500
    #     @info "Debugging" CACHE_SPAC.LATITUDE CACHE_SPAC.LONGITUDE GPP(CACHE_SPAC) PPAR(CACHE_SPAC);
    # end;

    return GPP(CACHE_SPAC)
);
