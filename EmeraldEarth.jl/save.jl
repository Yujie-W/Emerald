#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-14: add function to save result per time step
#
#######################################################################################################################################################################################################
"""

    save_simulations!(filename::String, states::Matrix{<:Union{Nothing,MultiLayerSPACState{FT}}}, doy::Number) where {FT}

Save the simulation results to netcdf file, given
- `filename` Path of the netcdf file
- `states` Matrix of states (including outputs to save)
- `doy` Day of year as a float

"""
function save_simulations!(filename::String, states::Matrix{<:Union{Nothing,MultiLayerSPACState{FT}}}, doy::Number) where {FT}
    # read results from matrix of states
    @inline get_value(state::Union{Nothing,MultiLayerSPACState}, fn::Symbol) = (
        return isnothing(state) ? NaN32 : Float32(getfield(state, fn));
    );
    _mat_gpp = get_value.(states, :gpp);
    _mat_sif₆₈₃ = get_value.(states, :tropomi_sif₆₈₃);
    _mat_sif₇₄₀ = get_value.(states, :tropomi_sif₇₄₀);

    # create file if not existing
    if !isfile(filename)
        create_nc!(filename, ["lon", "lat", "ind"], [size(states,1), size(states,2), Inf]);
        _res_lat = 180 / size(states,2);
        _res_lon = 360 / size(states,1);
        _3d_gpp = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_683 = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_740 = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_gpp[:,:,1] .= _mat_gpp;
        _3d_683[:,:,1] .= _mat_sif₆₈₃;
        _3d_740[:,:,1] .= _mat_sif₇₄₀;

        append_nc!(filename, "lat", collect(Float32, (-90+_res_lat):_res_lat:90), Dict{String,String}("about" => "latitude"), ["lat"]);
        append_nc!(filename, "lon", collect(Float32, (-180+_res_lon):_res_lon:180), Dict{String,String}("about" => "longitude"), ["lon"]);
        append_nc!(filename, "DOY", Float32[doy], Dict{String,String}("about" => "index of hour in a year"), ["ind"]);
        append_nc!(filename, "GPP", _3d_gpp, Dict{String,String}("about" => "GPP in [μmol m⁻² s⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "SIF₆₈₃", _3d_683, Dict{String,String}("about" => "TROPOMI SIF at 683 nm in [mW m⁻² nm⁻¹ sr⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "SIF₇₄₀", _3d_740, Dict{String,String}("about" => "TROPOMI SIF at 740 nm in [mW m⁻² nm⁻¹ sr⁻¹]"), ["lon", "lat", "ind"]);

        return nothing
    end;

    # grow the dataset
    grow_nc!(filename, "DOY", Float32(doy), true);
    grow_nc!(filename, "GPP", _mat_gpp, false);
    grow_nc!(filename, "SIF₆₈₃", _mat_sif₆₈₃, false);
    grow_nc!(filename, "SIF₇₄₀", _mat_sif₇₄₀, false);

    return nothing
end
