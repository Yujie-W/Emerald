#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-18: migrate from PkgUtility to JuliaUtility
#     2023-Feb-23: migrate from JuliaUtility to Emerald
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
