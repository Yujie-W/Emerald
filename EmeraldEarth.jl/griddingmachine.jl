#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: migrate from research repo to Emerald
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save gridded datasets from GriddingMachine

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LandDatasets{FT<:AbstractFloat}
    "Spatial resolution zoom factor, resolution is 1/gz °"
    gz::Int
    "Which year do the datasets apply (when applicable)"
    year::Int
    "GriddingMachine.jl tag for soil color class"
    tag_s_cc::String
    "GriddingMachine.jl tag for soil van Genuchten parameters"
    tag_s_α::String
    "GriddingMachine.jl tag for soil van Genuchten parameters"
    tag_s_n::String
    "GriddingMachine.jl tag for soil van Genuchten parameters"
    tag_s_Θr::String
    "GriddingMachine.jl tag for soil van Genuchten parameters"
    tag_s_Θs::String
    "GriddingMachine.jl tag for canopy height"
    tag_p_ch::String
    "GriddingMachine.jl tag for chlorophyll content"
    tag_p_chl::String
    "GriddingMachine.jl tag for clumping index"
    tag_p_ci::String
    "GriddingMachine.jl tag for leaf area index"
    tag_p_lai::String
    "GriddingMachine.jl tag for specific leaf area"
    tag_p_sla::String
    "GriddingMachine.jl tag for Vcmax"
    tag_p_vcm::String
    "GriddingMachine.jl tag for elevation"
    tag_t_ele::String
    "GriddingMachine.jl tag for land mask"
    tag_t_lm::String
    "GriddingMachine.jl tag for PFT"
    tag_t_pft::String

    # soil properties
    "Soil color class"
    s_cc::Array{FT} = regrid(read_LUT(query_collection(tag_s_cc))[1], gz);
    "Soil van Genuchten α"
    s_α::Array{FT} = regrid(read_LUT(query_collection(tag_s_α))[1], gz)
    "Soil van Genuchten n"
    s_n::Array{FT} = regrid(read_LUT(query_collection(tag_s_n))[1], gz)
    "Soil van Genuchten Θr"
    s_Θr::Array{FT} = regrid(read_LUT(query_collection(tag_s_Θr))[1], gz)
    "Soil van Genuchten Θs"
    s_Θs::Array{FT} = regrid(read_LUT(query_collection(tag_s_Θs))[1], gz)

    # plant properties
    "Plant canopy height"
    p_ch::Array{FT} = regrid(read_LUT(query_collection(tag_p_ch))[1], gz)
    "Plant chlorophyll content"
    p_chl::Array{FT} = regrid(read_LUT(query_collection(tag_p_chl))[1], gz)
    "Stand clumping index"
    p_ci::Array{FT} = regrid(read_LUT(query_collection(tag_p_ci))[1], gz)
    "Stand leaf area index"
    p_lai::Array{FT} = regrid(read_LUT(query_collection(tag_p_lai))[1], gz)
    "Plant leaf specific area"
    p_sla::Array{FT} = regrid(read_LUT(query_collection(tag_p_sla))[1], gz)
    "Plant maximum carboxylation rate"
    p_vcm::Array{FT} = regrid(read_LUT(query_collection(tag_p_vcm))[1], gz)

    # stand properties
    "Stand elevation"
    t_ele::Array{FT} = regrid(read_LUT(query_collection(tag_t_ele))[1], gz)
    "Stand land mask"
    t_lm::Array{FT} = regrid(read_LUT(query_collection(tag_t_lm))[1], gz)
    "Stand PFT percentages `[%]`"
    t_pft::Array{FT} = regrid(read_LUT(query_collection(tag_t_pft))[1], gz)
end

"""

    LandDatasets{FT}(gm_tag::String, year::Int) where {FT<:AbstractFloat}

Constructor of LandDatasets, given
- `gm_tag` Unique tag of GriddingMachine parameterization
- `year` year of simulations

"""
LandDatasets{FT}(gm_tag::String, year::Int) where {FT<:AbstractFloat} = (
    @assert gm_tag in ["gm1", "gm2"] "Parameterization tag $(gm_tag) is not supported!";

    if gm_tag == "gm1"
        return LandDatasets{FT}(
                    gz         = 1,
                    year       = year,
                    tag_s_cc   = "SC_2X_1Y_V1",
                    tag_s_α    = "SOIL_VGA_12X_1Y_V1",
                    tag_s_n    = "SOIL_VGN_12X_1Y_V1",
                    tag_s_Θr   = "SOIL_SWCR_12X_1Y_V1",
                    tag_s_Θs   = "SOIL_SWCS_12X_1Y_V1",
                    tag_p_ch   = "CH_20X_1Y_V1",
                    tag_p_chl  = "CHL_2X_7D_V1",
                    tag_p_ci   = "CI_2X_1Y_V1",
                    tag_p_lai  = "LAI_MODIS_2X_8D_$(year)_V1",
                    tag_p_sla  = "SLA_2X_1Y_V1",
                    tag_p_vcm  = "VCMAX_2X_1Y_V2",
                    tag_t_ele  = "ELEV_4X_1Y_V1",
                    tag_t_lm   = "LM_4X_1Y_V1",
                    tag_t_pft  = "PFT_2X_1Y_V1")
    elseif gm_tag == "gm2"
        return LandDatasets{FT}(
                    gz         = 1,
                    year       = year,
                    tag_s_cc   = "SC_2X_1Y_V1",
                    tag_s_α    = "SOIL_VGA_12X_1Y_V1",
                    tag_s_n    = "SOIL_VGN_12X_1Y_V1",
                    tag_s_Θr   = "SOIL_SWCR_12X_1Y_V1",
                    tag_s_Θs   = "SOIL_SWCS_12X_1Y_V1",
                    tag_p_ch   = "CH_20X_1Y_V1",
                    tag_p_chl  = "CHL_2X_7D_V1",
                    tag_p_ci   = "CI_2X_1M_V3",
                    tag_p_lai  = "LAI_MODIS_2X_8D_$(year)_V1",
                    tag_p_sla  = "SLA_2X_1Y_V1",
                    tag_p_vcm  = "VCMAX_2X_1Y_V2",
                    tag_t_ele  = "ELEV_4X_1Y_V1",
                    tag_t_lm   = "LM_4X_1Y_V1",
                    tag_t_pft  = "PFT_2X_1Y_V1")
    end;
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: migrate from research repo to Emerald
#
#######################################################################################################################################################################################################
"""

    extend_data!(data::Union{FT, Vector{FT}}) where {FT<:AbstractFloat}

Gap fill the data linearly, given
- `data` Input data

"""
function extend_data!(data::Union{FT, Vector{FT}}) where {FT<:AbstractFloat}
    if sum(.!isnan.(data)) in [0, length(data)]
        return nothing
    end;

    @inline find_last_number(vec_in::Vector{FT}, ind::Int) where {FT<:AbstractFloat} = (
        _x = ind;
        _y = vec_in[ind];
        for _i in ind:-1:1
            if !isnan(vec_in[_i])
                _x = _i;
                _y = vec_in[_i];
                break;
            end;
        end;

        return _x, _y
    );

    @inline find_next_number(vec_in::Vector{FT}, ind::Int) where {FT<:AbstractFloat} = (
        _x = ind;
        _y = vec_in[ind];
        for _i in ind:1:length(vec_in)
            if !isnan(vec_in[_i])
                _x = _i;
                _y = vec_in[_i];
                break;
            end;
        end;

        return _x, _y
    );

    @inline interpolate_data!(vec_in::Vector{FT}, ind::Int) where {FT<:AbstractFloat} = (
        if isnan(vec_in[ind])
            (_x1,_y1) = find_last_number(vec_in, ind);
            (_x2,_y2) = find_next_number(vec_in, ind);
            vec_in[ind] = ((ind - _x1) * _y2 + (_x2 - ind) * _y1) / (_x2 - _x1);
        end;

        return nothing
    );

    _data_3x = [data; data; data];
    interpolate_data!.([_data_3x], (length(data)+1):(length(data)*2));
    data .= _data_3x[(length(data)+1):(length(data)*2)];

    return nothing
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: migrate from research repo to Emerald
#
#######################################################################################################################################################################################################
"""

    query_grid_info_filename(gm_tag::String, year::Int)
    query_grid_info_filename(gm_tag::String, dts::LandDatasets)

Return the file path to grid information, given
- `gm_tag` Unique tag of GriddingMachine parameterization
- `year` Year of simulations
- `dts` `LandDatasets` type data to save

"""
function query_grid_info_filename end

query_grid_info_filename(gm_tag::String, year::Int) = "$(SETUP_FOLDER)/grid_info_$(gm_tag)_$(year).jld2";

query_grid_info_filename(gm_tag::String, dts::LandDatasets) = "$(SETUP_FOLDER)/grid_info_$(gm_tag)_$(dts.year).jld2";


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: migrate from research repo to Emerald
#
#######################################################################################################################################################################################################
"""

    prepare_grid_jld!(gm_tag::String, year::Int)
    prepare_grid_jld!(gm_tag::String, dts::LandDatasets)

Prepare the data used to run CliMA Land model, given
- `gm_tag` Unique tag of GriddingMachine parameterization (needs to be supported in LandDatasets constructor function)
- `year` Which year of data to read
- `dts` `LandDatasets` type data to read

"""
function prepare_grid_jld! end

prepare_grid_jld!(gm_tag::String, year::Int) = (
    # if file exists, do nothing
    _jld = query_grid_info_filename(gm_tag, year);
    if isfile(_jld)
        @info "File $(_jld) already exists!";
        return nothing
    end;

    # save the file if the file does not exist
    @info "Reading datasets for year $(year)...";
    _dts = LandDatasets{Float64}(gm_tag, year);

    prepare_grid_jld!(gm_tag, _dts);

    return nothing
);

prepare_grid_jld!(gm_tag::String, dts::LandDatasets) = (
    # if file exists, do nothing
    _jld = query_grid_info_filename(gm_tag, dts);
    _ccs = read_csv(artifact"2021_land_gpp" * "/CO2.csv");
    if isfile(_jld)
        @info "File $(_jld) already exists!";
        return nothing
    end;

    # combine lat and lon
    _res    = 1 / dts.gz;
    _lats   = collect(_res/2:_res:180) .- 90;
    _lons   = collect(_res/2:_res:360) .- 180;
    _ind_c3 = [2:14;16;17];
    _ind_c4 = [15];
    _dicts  = Dict{String,Any}[];
    for _lat in _lats, _lon in _lons
        _lat_ind = lat_ind(_lat; res=_res);
        _lon_ind = lon_ind(_lon; res=_res);
        _chls    = dts.p_chl[_lon_ind,_lat_ind,:];
        _ci      = dts.p_ci[_lon_ind,_lat_ind,:];
        _lais    = dts.p_lai[_lon_ind,_lat_ind,:];
        _lma     = 1 / dts.p_sla[_lon_ind,_lat_ind,1] / 10;
        _lmsk    = dts.t_lm[_lon_ind,_lat_ind,1];
        _pfts    = dts.t_pft[_lon_ind,_lat_ind,:];
        _scolor  = min(20, max(1, Int(floor(dts.s_cc[_lon_ind,_lat_ind,1]))));
        _s_α     = dts.s_α[_lon_ind,_lat_ind,:];
        _s_n     = dts.s_n[_lon_ind,_lat_ind,:];
        _s_Θr    = dts.s_Θr[_lon_ind,_lat_ind,:];
        _s_Θs    = dts.s_Θs[_lon_ind,_lat_ind,:];
        _vcmax   = dts.p_vcm[_lon_ind,_lat_ind,:];
        _zc      = dts.p_ch[_lon_ind,_lat_ind,1];

        # gap fill the data for seasonal trends
        extend_data!(_chls);
        extend_data!(_ci);
        extend_data!(_lais);
        extend_data!(_vcmax);

        # compute g1 for Medlyn model
        _g1_c3_medlyn = CLM5_PFTG[_ind_c3]' * _pfts[_ind_c3] / sum(_pfts[_ind_c3]);
        _g1_c4_medlyn = CLM5_PFTG[_ind_c4]' * _pfts[_ind_c4] / sum(_pfts[_ind_c4]);

        # CO2 concentration
        _co2 = _ccs.Mean[findfirst(_ccs.Year .== dts.year)];

        # filter the data by input data
        if (nanmax(_lais) > 0)      &&  # grid is vegetated
           (nanmax(_chls) > 0)      &&  # chlorophyll content is known
           (_lmsk > 0)              &&  # grid contains land
           (!isnan(_scolor))        &&  # soil color is known
           (all(.!isnan.(_vcmax)))  &&  # Vcmax is known
           (!isnan(_zc))            &&  # Canopy height is not NaN
           (!isnan(_lma))               # leaf mass per area is known
            _dict = Dict{String,Any}("canopy_height"      => _zc,
                                     "chlorophyll"        => _chls,
                                     "clumping_index"     => _ci,
                                     "co2_concentration"  => _co2,
                                     "g1_medlyn_c3"       => _g1_c3_medlyn,
                                     "g1_medlyn_c4"       => _g1_c4_medlyn,
                                     "land_mask"          => _lmsk,
                                     "latitude"           => _lat,
                                     "latitude_index"     => _lat_ind,
                                     "leaf_area_index"    => _lais,
                                     "leaf_mass_per_area" => _lma,
                                     "longitude"          => _lon,
                                     "longitude_index"    => _lon_ind,
                                     "pft_class"          => CLM5_PFTS,
                                     "pft_g1_medlyn"      => CLM5_PFTG,
                                     "pft_ratio"          => _pfts,
                                     "soil_color"         => _scolor,
                                     "soil_vg_n"          => _s_n,
                                     "soil_vg_α"          => _s_α,
                                     "soil_vg_Θr"         => _s_Θr,
                                     "soil_vg_Θs"         => _s_Θs,
                                     "spatial_resolution" => "$(dts.gz)X",
                                     "vcmax"              => _vcmax,
                                     "year"               => dts.year);
            push!(_dicts, _dict)
        end;
    end;

    # save the data to a new dictionary
    _clima_land_dict = Dict{String,Any}(
        "grid_info"             => _dicts,
        "griddingmachine_tags"  => String[dts.tag_s_cc, dts.tag_s_α, dts.tag_s_n, dts.tag_s_Θr, dts.tag_s_Θs,
                                          dts.tag_p_ch, dts.tag_p_chl, dts.tag_p_ci, dts.tag_p_lai, dts.tag_p_sla, dts.tag_p_vcm,
                                          dts.tag_t_ele, dts.tag_t_lm, dts.tag_t_pft],
    );
    save(_jld, _clima_land_dict);

    return nothing
);
