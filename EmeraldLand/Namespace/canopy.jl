#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-07: add CanopyOptics struct (will be a field for canopy structure)
#     2022-Jun-07: add more fields: fo, fs, po, ps, pso
#     2022-Jun-08: add more fields: ρ_dd, ρ_sd, σ_dob, σ_dof, σ_so, τ_dd, τ_sd, _ρ_dd, _ρ_sd, _τ_dd, _τ_sd
#     2022-Jun-09: rename variables to be more descriptive
#     2022-Jun-10: add more fields: p_sunlit, ϵ, ρ_lw, τ_lw, _mat_down, _mat_down, _ρ_lw, _τ_lw
#     2022-Jun-10: add more fields for sif calculations
#     2022-Jun-13: add more fields for sif calculations
#     2022-Jun-13: remove unnecessary cache variables
#     2022-Jun-15: rename to HyperspectralMLCanopyOpticalProperty
#     2022-Jul-19: use kwdef for the constructor
#     2022-Jul-19: add dimension control to struct
#     2023-May-19: use δlai per canopy layer
#     2023-Jun-16: remove fields DIM_*
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Verhoef LIDF algorithm

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralMLCanopyOpticalProperty{FT<:AbstractFloat}
    # Cache variables
    "Upwelling matrix for SIF excitation"
    _mat⁺::Matrix{FT}
    "Downwelling matrix for SIF excitation"
    _mat⁻::Matrix{FT}
    "Temporary cache used for matrix adding up purpose (DIM_INCL * DIM_AZI)"
    _tmp_mat_incl_azi_1::Matrix{FT}
    "Temporary cache used for vector operations (DIM_LAYER)"
    _tmp_vec_layer::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_1::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_2::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_3::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_4::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_5::Vector{FT}
    "Cache variable to store the SIF information"
    _tmp_vec_sif_6::Vector{FT}
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_1::Vector{FT}
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_2::Vector{FT}
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_3::Vector{FT}
end;

HyperspectralMLCanopyOpticalProperty(config::SPACConfiguration{FT}) where {FT} = (
    (; DIM_AZI, DIM_INCL, DIM_LAYER, DIM_SIF, DIM_SIFE, DIM_WL) = config;

    return HyperspectralMLCanopyOpticalProperty{FT}(
                _mat⁺               = zeros(FT, DIM_SIF, DIM_SIFE),
                _mat⁻               = zeros(FT, DIM_SIF, DIM_SIFE),
                _tmp_mat_incl_azi_1 = zeros(FT, DIM_INCL, DIM_AZI),
                _tmp_vec_layer      = zeros(FT, DIM_LAYER),
                _tmp_vec_sif_1      = zeros(FT, DIM_SIF),
                _tmp_vec_sif_2      = zeros(FT, DIM_SIF),
                _tmp_vec_sif_3      = zeros(FT, DIM_SIF),
                _tmp_vec_sif_4      = zeros(FT, DIM_SIF),
                _tmp_vec_sif_5      = zeros(FT, DIM_SIF),
                _tmp_vec_sif_6      = zeros(FT, DIM_SIF),
                _tmp_vec_sife_1     = zeros(FT, DIM_SIFE),
                _tmp_vec_sife_2     = zeros(FT, DIM_SIFE),
                _tmp_vec_sife_3     = zeros(FT, DIM_SIFE),
    )
);




#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-09: migrate CanopyRads as HyperspectralMLCanopyRadiationProfile
#     2022-Jun-09: add fields: albedo, apar_shaded, apar_sunlit, e_o, e_v, r_net
#     2022-Jun-10: add fields: r_net_sw, r_net_sw_shaded, r_net_sw_sunlit, r_lw, r_lw_down, r_lw_up, _r_emit_down, _r_emit_up
#     2022-Jun-10: add more fields for SIF
#     2022-Jun-13: add more fields for sif calculations
#     2022-Jun-15: rename to HyperspectralMLCanopyRadiationProfile
#     2022-Jun-27: move ppar_sunlit and ppar_shaded to Leaf
#     2022-Jul-19: add dimension control to struct
#     2022-Aug-30: rename sif_obs_ssoil to sif_obs_soil (typo fix, non-breaking)
#     2023-Jun-16: remove fields DIM_*
#     2023-Sep-11: add fields s_layer_down_chl and s_layer_up_chl to store the SIF at chloroplast level
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to store canopy radiation profiles

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralMLCanopyRadiationProfile{FT<:AbstractFloat}
    # Diagnostic variables
    "Albedo towards the viewing direction"
    albedo::Vector{FT}
    "Total radiation towards the viewing direction `[mW m⁻² nm⁻¹]`"
    e_o::Vector{FT}
    "Radiation towards the viewing direction per layer (including soil) `[mW m⁻² nm⁻¹]`"
    e_v::Matrix{FT}

    "Downwelling SIF for sunlit leaves at each wavelength for a layer"
    s_layer_down::Matrix{FT}
    "Downwelling SIF for sunlit leaves at each wavelength for a layer at chloroplast level"
    s_layer_down_chl::Matrix{FT}
    "Upwelling SIF for sunlit leaves at each wavelength for a layer"
    s_layer_up::Matrix{FT}
    "Upwelling SIF for sunlit leaves at each wavelength for a layer at chloroplast level"
    s_layer_up_chl::Matrix{FT}
    "Downwelling SIF"
    sif_down::Matrix{FT}
    "SIF at observer direction"
    sif_obs::Vector{FT}
    "SIF at observer direction from shaded APAR"
    sif_obs_shaded::Vector{FT}
    "SIF at observer direction from scattering"
    sif_obs_scatter::Vector{FT}
    "SIF at observer direction from soil reflection"
    sif_obs_soil::Vector{FT}
    "SIF at observer direction from sunlit APAR"
    sif_obs_sunlit::Vector{FT}
    "Upwelling SIF"
    sif_up::Matrix{FT}

    # Cache variables
    "Downwelling SIF for sunlit leaves at each wavelength"
    _s_emit_down::Matrix{FT}
    "Upwelling SIF for sunlit leaves at each wavelength"
    _s_emit_up::Matrix{FT}
    "Downwelling SIF for shaded leaves at each wavelength"
    _s_shaded_down::Vector{FT}
    "Upwelling SIF for shaded leaves at each wavelength"
    _s_shaded_up::Vector{FT}
    "Downwelling SIF for sunlit leaves at each wavelength"
    _s_sunlit_down::Vector{FT}
    "Upwelling SIF for sunlit leaves at each wavelength"
    _s_sunlit_up::Vector{FT}
    "Cache to compute SIF at observer direction from shaded APAR"
    _sif_obs_shaded::Matrix{FT}
    "Cache to compute SIF at observer direction from scattering"
    _sif_obs_scatter::Matrix{FT}
    "Cache to compute SIF at observer direction from sunlit APAR"
    _sif_obs_sunlit::Matrix{FT}
end;

HyperspectralMLCanopyRadiationProfile(config::SPACConfiguration{FT}) where {FT} = (
    (; DIM_AZI, DIM_INCL, DIM_LAYER, DIM_PAR, DIM_SIF, DIM_WL) = config;

    return HyperspectralMLCanopyRadiationProfile{FT}(
                albedo           = zeros(FT, DIM_WL),
                e_o              = zeros(FT, DIM_WL),
                e_v              = zeros(FT, DIM_WL, DIM_LAYER+1),
                s_layer_down     = zeros(FT, DIM_SIF, DIM_LAYER),
                s_layer_down_chl = zeros(FT, DIM_SIF, DIM_LAYER),
                s_layer_up       = zeros(FT, DIM_SIF, DIM_LAYER),
                s_layer_up_chl   = zeros(FT, DIM_SIF, DIM_LAYER),
                sif_down         = zeros(FT, DIM_SIF, DIM_LAYER+1),
                sif_obs          = zeros(FT, DIM_SIF),
                sif_obs_shaded   = zeros(FT, DIM_SIF),
                sif_obs_scatter  = zeros(FT, DIM_SIF),
                sif_obs_soil     = zeros(FT, DIM_SIF),
                sif_obs_sunlit   = zeros(FT, DIM_SIF),
                sif_up           = zeros(FT, DIM_SIF, DIM_LAYER+1),
                _s_emit_down     = zeros(FT, DIM_SIF, DIM_LAYER),
                _s_emit_up       = zeros(FT, DIM_SIF, DIM_LAYER+1),
                _s_shaded_down   = zeros(FT, DIM_SIF),
                _s_shaded_up     = zeros(FT, DIM_SIF),
                _s_sunlit_down   = zeros(FT, DIM_SIF),
                _s_sunlit_up     = zeros(FT, DIM_SIF),
                _sif_obs_shaded  = zeros(FT, DIM_SIF, DIM_LAYER),
                _sif_obs_scatter = zeros(FT, DIM_SIF, DIM_LAYER),
                _sif_obs_sunlit  = zeros(FT, DIM_SIF, DIM_LAYER),
    )
);




#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-02: migrate from CanopyLayers
#     2022-Jun-02: rename Canopy4RT to MultiLayerCanopy
#     2022-Jun-07: add cache variable _COS_Θ_INCL_AZI
#     2022-Jun-07: remove cache variable _vol_scatter
#     2022-Jun-09: add new field: APAR_CAR, RADIATION, WLSET
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jun-15: rename to HyperspectralMLCanopyOpticalProperty and HyperspectralMLCanopyRadiationProfile
#     2022-Jun-16: remove some cache variables
#     2022-Jul-22: remove field APAR_CAR
#     2022-Aug-30: add field LHA (moved from spac)
#     2023-May-22: add support to BetaLIDF
#     2023-Jun-15: compute _x_bnds based on lai (add case when lai = 0)
#     2023-Jun-16: move field WLSET to SPACConfiguration
#     2023-Jun-16: remove fields DIM_*
#     2023-Jun-20: move LHA to SPACConfiguration
#     2023-Jun-20: move fields Θ_AZI, Θ_INCL, Θ_INCL_BNDS, _COS_Θ_INCL_AZI to SPACConfiguration
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save multiple layer hyperspectral canopy parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MultiLayerCanopy{FT<:AbstractFloat}
    # sun-sensor geometry related structs
    "Sensor geometry information"
    sensor_geometry::SensorGeometry{FT}
    "Sun geometry information"
    sun_geometry::SunGeometry{FT}
    "Canopy structure"
    structure::CanopyStructure{FT}




    # General model information

    # Embedded structures
    "Canopy optical properties"
    OPTICS::HyperspectralMLCanopyOpticalProperty{FT}
    "Canopy radiation profiles"
    RADIATION::HyperspectralMLCanopyRadiationProfile{FT}

    # Geometry information
end;

MultiLayerCanopy(config::SPACConfiguration{FT}) where {FT} = (
    return MultiLayerCanopy{FT}(
                sensor_geometry = SensorGeometry(config),
                sun_geometry    = SunGeometry(config),
                structure       = CanopyStructure(config),
                OPTICS          = HyperspectralMLCanopyOpticalProperty(config),
                RADIATION       = HyperspectralMLCanopyRadiationProfile(config),
    )
);
