#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-14: add function to save result per time step
#
#######################################################################################################################################################################################################
"""

    save_simulations!(filename::String, states::Matrix{Union{Nothing,MultiLayerSPACState{FT}}}, doy::Number) where {FT<:AbstractFloat}

Save the simulation results to netcdf file, given
- `filename` Path of the netcdf file
- `states` Matrix of states (including outputs to save)
- `doy` Day of year as a float

"""
function save_simulations!(filename::String, states::Matrix{Union{Nothing,MultiLayerSPACState{FT}}}, doy::Number) where {FT<:AbstractFloat}
    # read results from matrix of states
    @inline get_value(state::Union{Nothing,MultiLayerSPACState}, fn::Symbol) = (
        return isnothing(state) ? NaN32 : Float32(getfield(state, fn));
    );
    _mat_gpp = get_value.(states, :gpp);
    _mat_evi = get_value.(states, :modis_evi);
    _mat_ndvi = get_value.(states, :modis_ndvi);
    _mat_nirv = get_value.(states, :modis_nirv);
    _mat_sif₆₈₃ = get_value.(states, :tropomi_sif₆₈₃);
    _mat_sif₇₄₀ = get_value.(states, :tropomi_sif₇₄₀);
    _mat_sif₇₅₉ = get_value.(states, :oco_sif₇₅₉);
    _mat_sif₇₇₀ = get_value.(states, :oco_sif₇₇₀);

    # create file if not existing
    if !isfile(filename)
        create_nc!(filename, ["lon", "lat", "ind"], [size(states,1), size(states,2), Inf]);
        _res_lat = 180 / size(states,2);
        _res_lon = 360 / size(states,1);
        _3d_gpp = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_evi = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_ndvi = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_nirv = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_683 = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_740 = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_759 = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_770 = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_gpp[:,:,1] .= _mat_gpp;
        _3d_evi[:,:,1] .= _mat_evi;
        _3d_ndvi[:,:,1] .= _mat_ndvi;
        _3d_nirv[:,:,1] .= _mat_nirv;
        _3d_683[:,:,1] .= _mat_sif₆₈₃;
        _3d_740[:,:,1] .= _mat_sif₇₄₀;
        _3d_759[:,:,1] .= _mat_sif₇₅₉;
        _3d_770[:,:,1] .= _mat_sif₇₇₀;

        append_nc!(filename, "lat", collect(Float32, (-90+_res_lat):_res_lat:90), Dict{String,String}("about" => "latitude"), ["lat"]);
        append_nc!(filename, "lon", collect(Float32, (-180+_res_lon):_res_lon:180), Dict{String,String}("about" => "longitude"), ["lon"]);
        append_nc!(filename, "DOY", Float32[doy], Dict{String,String}("about" => "index of hour in a year"), ["ind"]);
        append_nc!(filename, "GPP", _3d_gpp, Dict{String,String}("about" => "GPP in [μmol m⁻² s⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "EVI", _3d_evi, Dict{String,String}("about" => "MODIS EVI"), ["lon", "lat", "ind"]);
        append_nc!(filename, "NDVI", _3d_ndvi, Dict{String,String}("about" => "MODIS NDVI"), ["lon", "lat", "ind"]);
        append_nc!(filename, "NIRv", _3d_nirv, Dict{String,String}("about" => "MODIS NIRv"), ["lon", "lat", "ind"]);
        append_nc!(filename, "SIF₆₈₃", _3d_683, Dict{String,String}("about" => "TROPOMI SIF at 683 nm in [mW m⁻² nm⁻¹ sr⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "SIF₇₄₀", _3d_740, Dict{String,String}("about" => "TROPOMI SIF at 740 nm in [mW m⁻² nm⁻¹ sr⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "SIF₇₅₉", _3d_759, Dict{String,String}("about" => "OCO SIF at 759 nm in [mW m⁻² nm⁻¹ sr⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "SIF₇₇₀", _3d_770, Dict{String,String}("about" => "OCO SIF at 770 nm in [mW m⁻² nm⁻¹ sr⁻¹]"), ["lon", "lat", "ind"]);

        return nothing
    end;

    # grow the dataset
    grow_nc!(filename, "DOY", Float32(doy), true);
    grow_nc!(filename, "GPP", _mat_gpp, false);
    grow_nc!(filename, "EVI", _mat_evi, false);
    grow_nc!(filename, "NDVI", _mat_ndvi, false);
    grow_nc!(filename, "NIRv", _mat_nirv, false);
    grow_nc!(filename, "SIF₆₈₃", _mat_sif₆₈₃, false);
    grow_nc!(filename, "SIF₇₄₀", _mat_sif₇₄₀, false);
    grow_nc!(filename, "SIF₇₅₉", _mat_sif₇₅₉, false);
    grow_nc!(filename, "SIF₇₇₀", _mat_sif₇₇₀, false);

    return nothing
end
