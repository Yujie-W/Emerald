#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to run the SPAC
#     2023-Mar-13: add state mat as an input
#     2023-Apr-13: add CACHE_CONFIG in function
#
#######################################################################################################################################################################################################
"""

    simulation!(gm_mat::Matrix{Union{Nothing,Dict{String,Any}}}, wd_mat::Matrix{Union{Nothing,Dict{String,Any}}}, state_mat::Matrix{Union{Nothing,MultiLayerSPACState{FT}}}) where {FT<:AbstractFloat}

Run simulations on SPAC, given
- `gm_mat` Matrix of GriddingMachine inputs
- `wd_mat` Matrix of weather drivers
- `state_mat` Matrix of state variable struct

---
# Example
```julia
using Emerald;

FT = Float64;
@time EmeraldEarth.add_threads!(20, FT);

@time dts = EmeraldEarth.LandDatasets{FT}("gm2", 2020);
@time mat = EmeraldEarth.gm_grids(dts);
@time sts = Matrix{Union{Nothing,EmeraldCore.Namespace.MultiLayerSPACState{FT}}}(nothing, size(dts.t_lm));
@time wdr = EmeraldEarth.ERA5SingleLevelsDriver();

@time wds = EmeraldEarth.weather_drivers(dts, wdr);
@time wdx = EmeraldEarth.wd_grids(dts, wds, 1);

# for debugging use
@time wdx = EmeraldEarth.wd_grids(dts, wdr, 1);
@time sts = EmeraldEarth.simulation!(mat, wdx, sts);
@time EmeraldEarth.save_simulations!("test.nc", sts, 1);
@time wdx = EmeraldEarth.wd_grids(dts, wdr, 2);
@time sts = EmeraldEarth.simulation!(mat, wdx, sts);
@time EmeraldEarth.save_simulations!("test.nc", sts, 2);

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

simulation!(gm_mat::Matrix{Union{Nothing,Dict{String,Any}}}, wd_mat::Matrix{Union{Nothing,Dict{String,Any}}}, state_mat::Matrix{Union{Nothing,MultiLayerSPACState{FT}}}) where {FT<:AbstractFloat} = (
    @tinfo "Running the global simulations in multiple threads...";
    _states = @showprogress pmap(simulation!, gm_mat, wd_mat, state_mat);

    return _states
);

simulation!(gm_params::Nothing, wd_params::Nothing, state::Nothing) = nothing;

simulation!(gm_params::Dict{String,Any}, wd_params::Dict{String,Any}, state::Union{Nothing,MultiLayerSPACState{FT}}) where {FT<:AbstractFloat} = (
    synchronize_cache!(gm_params, wd_params, state);
    for _ in 1:10
        soil_plant_air_continuum!(CACHE_SPAC, CACHE_CONFIG, 360; p_on = false, t_on = false, θ_on = false);
    end;
    spac_state!(CACHE_SPAC, CACHE_STATE, CACHE_CONFIG);

    # if wd_params["RAD_DIR"] > 500
    #     @info "Debugging" CACHE_SPAC.LATITUDE CACHE_SPAC.LONGITUDE GPP(CACHE_SPAC) PPAR(CACHE_SPAC);
    # end;

    return CACHE_STATE
);
