#######################################################################################################################################################################################################
#
# Changes to this struct
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

    # masks
    "Mask for bare soil"
    mask_soil::Matrix{Bool} = zeros(Bool, size(t_lm))
    "Mask for SPAC"
    mask_spac::Matrix{Bool} = zeros(Bool, size(t_lm))
end

"""

    LandDatasets{FT}(gm_tag::String, year::Int) where {FT<:AbstractFloat}

Constructor of LandDatasets, given
- `gm_tag` Unique tag of GriddingMachine parameterization
- `year` year of simulations

"""
LandDatasets{FT}(gm_tag::String, year::Int) where {FT<:AbstractFloat} = (
    @assert gm_tag in ["gm1", "gm2"] "Parameterization tag $(gm_tag) is not supported!";

    @tinfo "Querying data from GriddingMachine...";
    if gm_tag == "gm1"
        _dts = LandDatasets{FT}(
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
    else # gm_tag == "gm2"
        _dts = LandDatasets{FT}(
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

    @tinfo "Gap-filling data from GriddingMachine...";
    extend_data!(_dts);

    return _dts
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: migrate from research repo to Emerald
#     2023-Mar-11: add method to extend the LandDatasets
#
#######################################################################################################################################################################################################
"""

    extend_data!(dts::LandDatasets{FT}) where {FT<:AbstractFloat}
    extend_data!(data::Union{FT, Vector{FT}}) where {FT<:AbstractFloat}

Gap fill the data linearly, given
- `dts` LandDatasets struct
- `data` Input data

"""
function extend_data! end

extend_data!(dts::LandDatasets{FT}) where {FT<:AbstractFloat} = (
    # determine where to fill based on land mask and lai
    for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
        if (dts.t_lm[_ilon,_ilat] > 0) && (nanmax(dts.p_lai[_ilon,_ilat,:]) > 0)
            dts.mask_spac[_ilon,_ilat] = true;
            _mask_lai = isnan.(dts.p_lai[_ilon,_ilat,:]);
            dts.p_lai[_ilon,_ilat,_mask_lai] .= 0;
        elseif (dts.t_lm[_ilon,_ilat] > 0)
            dts.mask_soil[_ilon,_ilat] = true;
            _mask_lai = isnan.(dts.p_lai[_ilon,_ilat,:]);
            dts.p_lai[_ilon,_ilat,_mask_lai] .= 0;
        end;
    end;

    # iterate the fieldnames
    for _field in fieldnames(typeof(dts))
        if !(_field in [:p_lai, :t_ele, :t_lm, :t_pft, :mask_soil, :mask_spac])
            _data = getfield(dts, _field);
            if _data isa Array
                # extend the data first based on interpolations
                for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
                    _tmp = _data[_ilon,_ilat,:];
                    extend_data!(_tmp);
                    _data[_ilon,_ilat,:] .= _tmp;
                end;

                # fill the NaNs with nanmean of the rest
                _mask_mean = dts.mask_spac .&& isnan.(_data);
                _data[_mask_mean] .= nanmean(_data);
            end;
        end;
    end;

    return nothing
);

extend_data!(data::Union{FT, Vector{FT}}) where {FT<:AbstractFloat} = (
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
);


# CLM5 settings
CLM5_PFTG = [0, 2.35, 2.35, 2.35, 4.12, 4.12, 4.45, 4.45, 4.45, 4.7, 4.7, 4.7, 2.22, 5.25, 1.62, 5.79, 5.79] .* sqrt(1000);
CLM5_PFTS = ["not_vegetated",
             "needleleaf_evergreen_temperate",
             "needleleaf_evergreen_boreal",
             "needleleaf_deciduous_boreal",
             "broadleaf_evergreen_tropical",
             "broadleaf_evergreen_temperate",
             "broadleaf_deciduous_tropical",
             "broadleaf_deciduous_temperate",
             "broadleaf_deciduous_boreal",
             "evergreen_shrub",
             "deciduous_temperate_shrub",
             "deciduous_boreal_shrub",
             "c3_arctic_grass",
             "c3_non-arctic_grass",
             "c4_grass",
             "c3_crop",
             "c3_irrigated"];


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: migrate from research repo to Emerald
#     2023-Jun-15: add non-vegetated land in global simulations
#
#######################################################################################################################################################################################################
"""

    gm_grids(dts::LandDatasets{FT}) where {FT<:AbstractFloat}

Prepare a matrix of GriddingMachine data to feed SPAC, given
- `dts` `LandDatasets` type data struct

"""
function gm_grids(dts::LandDatasets{FT}) where {FT<:AbstractFloat}
    # read some general data
    _ind_c3 = [2:14;16;17];
    _ccs = read_csv("$(@__DIR__)/../data/CO2-1Y.csv");
    _co2 = _ccs.MEAN[findfirst(_ccs.YEAR .== dts.year)];

    # create a matrix of GriddingMachine data
    # TODO: add a step to verify the input datasets
    @tinfo "Preparing a matrix of GriddingMachine data to work on...";
    _mat_gm = Matrix{Union{Nothing,Dict{String,Any}}}(nothing, size(dts.t_lm));
    for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
        if dts.mask_spac[_ilon,_ilat]
            _pfts = dts.t_pft[_ilon,_ilat,:];
            _g = CLM5_PFTG[_ind_c3]' * _pfts[_ind_c3] / sum(_pfts[_ind_c3]);
            _g1 = isnan(_g) ? nanmean(CLM5_PFTG[_ind_c3]) : _g;
            _mat_gm[_ilon,_ilat] = Dict{String,Any}(
                        "CANOPY_HEIGHT" => dts.p_ch[_ilon,_ilat],
                        "CHLOROPHYLL"   => dts.p_chl[_ilon,_ilat,:],
                        "CLUMPING"      => dts.p_ci[_ilon,_ilat,:],
                        "CO2"           => _co2,
                        "ELEVATION"     => dts.t_ele[_ilon,_ilat],
                        "FT"            => FT,
                        "LAI"           => dts.p_lai[_ilon,_ilat,:],
                        "LATITUDE"      => (_ilat - 0.5) * 180 / size(dts.t_lm,2) - 90,
                        "LMA"           => 1 / dts.p_sla[_ilon,_ilat] / 10,
                        "LONGITUDE"     => (_ilon - 0.5) * 360 / size(dts.t_lm,1) - 180,
                        "MEDLYN_G1"     => _g1,
                        "SOIL_COLOR"    => min(20, max(1, Int(floor(dts.s_cc[_ilon,_ilat])))),
                        "SOIL_N"        => dts.s_n[_ilon,_ilat,:],
                        "SOIL_α"        => dts.s_α[_ilon,_ilat,:],
                        "SOIL_ΘR"       => dts.s_Θr[_ilon,_ilat,:],
                        "SOIL_ΘS"       => dts.s_Θs[_ilon,_ilat,:],
                        "VCMAX25"       => dts.p_vcm[_ilon,_ilat,:],
                        "YEAR"          => dts.year,
            );
        elseif dts.mask_soil[_ilon,_ilat]
            _mat_gm[_ilon,_ilat] = Dict{String,Any}(
                        "CANOPY_HEIGHT" => eps(FT),
                        "CHLOROPHYLL"   => [10.0],
                        "CLUMPING"      => [1.0],
                        "CO2"           => _co2,
                        "ELEVATION"     => dts.t_ele[_ilon,_ilat],
                        "FT"            => FT,
                        "LAI"           => [0.0],
                        "LATITUDE"      => (_ilat - 0.5) * 180 / size(dts.t_lm,2) - 90,
                        "LMA"           => 0.012,
                        "LONGITUDE"     => (_ilon - 0.5) * 360 / size(dts.t_lm,1) - 180,
                        "MEDLYN_G1"     => nanmean(CLM5_PFTG[_ind_c3]),
                        "SOIL_COLOR"    => min(20, max(1, Int(floor(dts.s_cc[_ilon,_ilat])))),
                        "SOIL_N"        => dts.s_n[_ilon,_ilat,:],  # these are nan for some sites, fix this with a better data source
                        "SOIL_α"        => dts.s_α[_ilon,_ilat,:],  # these are nan for some sites, fix this with a better data source
                        "SOIL_ΘR"       => dts.s_Θr[_ilon,_ilat,:], # these are nan for some sites, fix this with a better data source
                        "SOIL_ΘS"       => dts.s_Θs[_ilon,_ilat,:], # these are nan for some sites, fix this with a better data source
                        "VCMAX25"       => [10.0],
                        "YEAR"          => dts.year,
            );
        end;
    end;

    return _mat_gm
end
