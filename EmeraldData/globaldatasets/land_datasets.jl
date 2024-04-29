#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-22: separate the labels from the datasets
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save gridded datasets from GriddingMachine (labels only)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LandDatasetLabels
    "GriddingMachine datasets version"
    gm_tag::String
    "Which year do the datasets apply (when applicable)"
    year::Int
    "Spatial resolution zoom factor, resolution is 1/nx °"
    nx::Int
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
end;

"""

    LandDatasetLabels(gm_tag::String, year::Int)

Constructor of LandDatasetLabels, given
- `gm_tag` version of GriddingMachine datasets collection
- `year` year of the datasets

"""
LandDatasetLabels(gm_tag::String, year::Int) = (
    @assert gm_tag in ["gm1", "gm2", "gm3", "gm4"] "Parameterization tag $(gm_tag) is not supported!";

    if gm_tag == "gm1"
        dtl = LandDatasetLabels(
                    gm_tag    = gm_tag,
                    year      = year,
                    nx        = 1,
                    tag_s_cc  = "SC_2X_1Y_V1",
                    tag_s_α   = "SOIL_VGA_12X_1Y_V1",
                    tag_s_n   = "SOIL_VGN_12X_1Y_V1",
                    tag_s_Θr  = "SOIL_SWCR_12X_1Y_V1",
                    tag_s_Θs  = "SOIL_SWCS_12X_1Y_V1",
                    tag_p_ch  = "CH_20X_1Y_V1",
                    tag_p_chl = "CHL_2X_7D_V1",
                    tag_p_ci  = "CI_2X_1Y_V1",
                    tag_p_lai = "LAI_MODIS_2X_8D_$(year)_V1",
                    tag_p_sla = "SLA_2X_1Y_V1",
                    tag_p_vcm = "VCMAX_2X_1Y_V2",
                    tag_t_ele = "ELEV_4X_1Y_V1",
                    tag_t_lm  = "LM_4X_1Y_V1",
                    tag_t_pft = "PFT_2X_1Y_V1")
    elseif gm_tag == "gm2"
        dtl = LandDatasetLabels(
                    gm_tag    = gm_tag,
                    year      = year,
                    nx        = 1,
                    tag_s_cc  = "SC_2X_1Y_V1",
                    tag_s_α   = "SOIL_VGA_12X_1Y_V1",
                    tag_s_n   = "SOIL_VGN_12X_1Y_V1",
                    tag_s_Θr  = "SOIL_SWCR_12X_1Y_V1",
                    tag_s_Θs  = "SOIL_SWCS_12X_1Y_V1",
                    tag_p_ch  = "CH_20X_1Y_V1",
                    tag_p_chl = "CHL_2X_7D_V1",
                    tag_p_ci  = "CI_2X_1M_V3",
                    tag_p_lai = "LAI_MODIS_2X_8D_$(year)_V1",
                    tag_p_sla = "SLA_2X_1Y_V1",
                    tag_p_vcm = "VCMAX_2X_1Y_V2",
                    tag_t_ele = "ELEV_4X_1Y_V1",
                    tag_t_lm  = "LM_4X_1Y_V1",
                    tag_t_pft = "PFT_2X_1Y_V1")
    elseif gm_tag == "gm3"
        dtl = LandDatasetLabels(
                    gm_tag    = gm_tag,
                    year      = year,
                    nx        = 1,
                    tag_s_cc  = "SC_2X_1Y_V1",
                    tag_s_α   = "SOIL_VGA_12X_1Y_V1",
                    tag_s_n   = "SOIL_VGN_12X_1Y_V1",
                    tag_s_Θr  = "SOIL_SWCR_12X_1Y_V1",
                    tag_s_Θs  = "SOIL_SWCS_12X_1Y_V1",
                    tag_p_ch  = "CH_20X_1Y_V1",
                    tag_p_chl = "CHL_2X_7D_V1",
                    tag_p_ci  = "CI_2X_1M_V3",
                    tag_p_lai = "LAI_MODIS_2X_8D_$(year)_V1",
                    tag_p_sla = "SLA_2X_1Y_V1",
                    tag_p_vcm = "VCMAX_2X_1Y_V2_60%",
                    tag_t_ele = "ELEV_4X_1Y_V1",
                    tag_t_lm  = "LM_4X_1Y_V1",
                    tag_t_pft = "PFT_2X_1Y_V1")
    elseif gm_tag == "gm4"
        dtl = LandDatasetLabels(
                    gm_tag    = gm_tag,
                    year      = year,
                    nx        = 1,
                    tag_s_cc  = "SC_2X_1Y_V1",
                    tag_s_α   = "SOIL_VGA_12X_1Y_V1",
                    tag_s_n   = "SOIL_VGN_12X_1Y_V1",
                    tag_s_Θr  = "SOIL_SWCR_12X_1Y_V1",
                    tag_s_Θs  = "SOIL_SWCS_12X_1Y_V1",
                    tag_p_ch  = "CH_20X_1Y_V1",
                    tag_p_chl = "CHL_2X_7D_V1",
                    tag_p_ci  = "CI_2X_1M_V3",
                    tag_p_lai = "LAI_MODIS_2X_8D_$(year)_V1",
                    tag_p_sla = "SLA_2X_1Y_V1",
                    tag_p_vcm = "VCMAX_CESM_1X_1M_V3",
                    tag_t_ele = "ELEV_4X_1Y_V1",
                    tag_t_lm  = "LM_4X_1Y_V1",
                    tag_t_pft = "PFT_2X_1Y_V1")
    else
        error("Tag $(gm_tag) is not supported!");
    end;

    return dtl
);


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
    # labels for the datasets
    "Land dataset labels"
    LABELS::LandDatasetLabels

    # soil properties
    "Soil color class"
    s_cc::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_s_cc))[1], LABELS.nx);
    "Soil van Genuchten α"
    s_α::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_s_α))[1], LABELS.nx)
    "Soil van Genuchten n"
    s_n::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_s_n))[1], LABELS.nx)
    "Soil van Genuchten Θr"
    s_Θr::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_s_Θr))[1], LABELS.nx)
    "Soil van Genuchten Θs"
    s_Θs::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_s_Θs))[1], LABELS.nx)

    # plant properties
    "Plant canopy height"
    p_ch::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_p_ch))[1], LABELS.nx)
    "Plant chlorophyll content"
    p_chl::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_p_chl))[1], LABELS.nx)
    "Stand clumping index"
    p_ci::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_p_ci))[1], LABELS.nx)
    "Stand leaf area index"
    p_lai::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_p_lai))[1], LABELS.nx)
    "Plant leaf specific area"
    p_sla::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_p_sla))[1], LABELS.nx)
    "Plant maximum carboxylation rate"
    p_vcm::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_p_vcm))[1], LABELS.nx)

    # stand properties
    "Stand elevation"
    t_ele::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_t_ele))[1], LABELS.nx)
    "Stand land mask"
    t_lm::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_t_lm))[1], LABELS.nx)
    "Stand PFT percentages `[%]`"
    t_pft::Array{FT} = regrid(read_LUT(query_collection(LABELS.tag_t_pft))[1], LABELS.nx)

    # masks
    "Mask for bare soil"
    mask_soil::Matrix{Bool} = zeros(Bool, size(t_lm))
    "Mask for SPAC"
    mask_spac::Matrix{Bool} = zeros(Bool, size(t_lm))
end;

"""

    LandDatasets{FT}(gm_tag::String, year::Int) where {FT}

Constructor of LandDatasets, given
- `gm_tag` Unique tag of GriddingMachine parameterization
- `year` year of simulations

"""
LandDatasets{FT}(gm_tag::String, year::Int) where {FT} = (
    dtl = LandDatasetLabels(gm_tag, year);
    @tinfo "Querying data from GriddingMachine...";
    if gm_tag in ["gm1", "gm2", "gm4"]
        dts = LandDatasets{FT}(LABELS = dtl);
    elseif gm_tag == "gm3"
        dts = LandDatasets{FT}(LABELS = dtl, p_vcm = regrid(read_LUT(query_collection("VCMAX_2X_1Y_V2"))[1], 1) .* 0.6);
    else
        error("Tag $(gm_tag) is not supported!");
    end;

    @tinfo "Gap-filling data from GriddingMachine...";
    extend_data!(dts);

    return dts
);
