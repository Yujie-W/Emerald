#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-14: add function to save result per time step
#     2023-Jun-15: save more data in the netcdf file
#     2023-Jun-15: add option to display process information
#
#######################################################################################################################################################################################################
"""

    save_simulations!(filename::String, states::Matrix{Union{Nothing,MultiLayerSPACState{FT}}}, doy::Number; displaying::Bool = true) where {FT}

Save the simulation results to netcdf file, given
- `filename` Path of the netcdf file
- `states` Matrix of states (including outputs to save)
- `doy` Day of year as a float
- `displaying` Whether to display information regarding process

"""
function save_simulations!(filename::String, states::Matrix{Union{Nothing,MultiLayerSPACState{FT}}}, doy::Number; displaying::Bool = true) where {FT}
    if displaying
        @tinfo "Saving the simulation results to netcdf file...";
    end;

    # read results from matrix of states
    @inline get_value(state::Union{Nothing,MultiLayerSPACState}, fn::Symbol) = (
        return isnothing(state) ? NaN32 : Float32(getfield(state, fn));
    );
    _mat_beta = get_value.(states, :beta);
    _mat_csif = get_value.(states, :csif);
    _mat_etr = get_value.(states, :etr);
    _mat_gpp = get_value.(states, :gpp);
    _mat_evi = get_value.(states, :modis_evi);
    _mat_ndvi = get_value.(states, :modis_ndvi);
    _mat_nirv = get_value.(states, :modis_nirv);
    _mat_par = get_value.(states, :par);
    _mat_ppar = get_value.(states, :ppar);
    _mat_sif₆₈₃ = get_value.(states, :tropomi_sif₆₈₃);
    _mat_sif₇₄₀ = get_value.(states, :tropomi_sif₇₄₀);
    _mat_sif₇₅₉ = get_value.(states, :oco_sif₇₅₉);
    _mat_sif₇₇₀ = get_value.(states, :oco_sif₇₇₀);
    _mat_tran = get_value.(states, :transpiration);

    # create file if not existing
    if !isfile(filename)
        create_nc!(filename, ["lon", "lat", "ind"], [size(states,1), size(states,2), Inf]);
        _res_lat = 180 / size(states,2);
        _res_lon = 360 / size(states,1);
        _3d_beta = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_csif = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_etr = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_gpp = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_evi = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_ndvi = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_nirv = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_par = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_ppar = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_s683 = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_s740 = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_s759 = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_s770 = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_tran = ones(Float32, size(_mat_gpp,1), size(_mat_gpp,2), 1);
        _3d_beta[:,:,1] .= _mat_beta;
        _3d_csif[:,:,1] .= _mat_csif;
        _3d_etr[:,:,1] .= _mat_etr;
        _3d_gpp[:,:,1] .= _mat_gpp;
        _3d_evi[:,:,1] .= _mat_evi;
        _3d_ndvi[:,:,1] .= _mat_ndvi;
        _3d_nirv[:,:,1] .= _mat_nirv;
        _3d_par[:,:,1] .= _mat_par;
        _3d_ppar[:,:,1] .= _mat_ppar;
        _3d_s683[:,:,1] .= _mat_sif₆₈₃;
        _3d_s740[:,:,1] .= _mat_sif₇₄₀;
        _3d_s759[:,:,1] .= _mat_sif₇₅₉;
        _3d_s770[:,:,1] .= _mat_sif₇₇₀;
        _3d_tran[:,:,1] .= _mat_tran;

        append_nc!(filename, "lat", collect(Float32, (-90+_res_lat):_res_lat:90), Dict{String,String}("about" => "latitude"), ["lat"]);
        append_nc!(filename, "lon", collect(Float32, (-180+_res_lon):_res_lon:180), Dict{String,String}("about" => "longitude"), ["lon"]);
        append_nc!(filename, "DOY", Float32[doy], Dict{String,String}("about" => "index of hour in a year"), ["ind"]);
        append_nc!(filename, "BETA", _3d_beta, Dict{String,String}("about" => "Beta factor"), ["lon", "lat", "ind"]);
        append_nc!(filename, "CSIF", _3d_csif, Dict{String,String}("about" => "Total Chlorophyll SIF photons in [μmol m⁻² s⁻¹] (before reabsorption)"), ["lon", "lat", "ind"]);
        append_nc!(filename, "ETR", _3d_etr, Dict{String,String}("about" => "Total electron transport in [μmol m⁻² s⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "GPP", _3d_gpp, Dict{String,String}("about" => "GPP in [μmol m⁻² s⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "EVI", _3d_evi, Dict{String,String}("about" => "MODIS EVI"), ["lon", "lat", "ind"]);
        append_nc!(filename, "NDVI", _3d_ndvi, Dict{String,String}("about" => "MODIS NDVI"), ["lon", "lat", "ind"]);
        append_nc!(filename, "NIRv", _3d_nirv, Dict{String,String}("about" => "MODIS NIRv"), ["lon", "lat", "ind"]);
        append_nc!(filename, "PAR", _3d_par, Dict{String,String}("about" => "Incoming PAR in [μmol m⁻² s⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "PPAR", _3d_ppar, Dict{String,String}("about" => "Photosynthesis PAR in [μmol m⁻² s⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "SIF₆₈₃", _3d_s683, Dict{String,String}("about" => "TROPOMI SIF at 683 nm in [mW m⁻² nm⁻¹ sr⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "SIF₇₄₀", _3d_s740, Dict{String,String}("about" => "TROPOMI SIF at 740 nm in [mW m⁻² nm⁻¹ sr⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "SIF₇₅₉", _3d_s759, Dict{String,String}("about" => "OCO SIF at 759 nm in [mW m⁻² nm⁻¹ sr⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "SIF₇₇₀", _3d_s770, Dict{String,String}("about" => "OCO SIF at 770 nm in [mW m⁻² nm⁻¹ sr⁻¹]"), ["lon", "lat", "ind"]);
        append_nc!(filename, "TRANSPIRATION", _3d_tran, Dict{String,String}("about" => "Transpiration rate in [mol m⁻² s⁻¹]"), ["lon", "lat", "ind"]);

        return nothing
    end;

    # grow the dataset
    grow_nc!(filename, "DOY", Float32(doy), true);
    grow_nc!(filename, "BETA", _mat_beta, false);
    grow_nc!(filename, "CSIF", _mat_csif, false);
    grow_nc!(filename, "ETR", _mat_etr, false);
    grow_nc!(filename, "GPP", _mat_gpp, false);
    grow_nc!(filename, "EVI", _mat_evi, false);
    grow_nc!(filename, "NDVI", _mat_ndvi, false);
    grow_nc!(filename, "NIRv", _mat_nirv, false);
    grow_nc!(filename, "PAR", _mat_par, false);
    grow_nc!(filename, "PPAR", _mat_ppar, false);
    grow_nc!(filename, "SIF₆₈₃", _mat_sif₆₈₃, false);
    grow_nc!(filename, "SIF₇₄₀", _mat_sif₇₄₀, false);
    grow_nc!(filename, "SIF₇₅₉", _mat_sif₇₅₉, false);
    grow_nc!(filename, "SIF₇₇₀", _mat_sif₇₇₀, false);
    grow_nc!(filename, "TRANSPIRATION", _mat_tran, false);

    return nothing
end;
