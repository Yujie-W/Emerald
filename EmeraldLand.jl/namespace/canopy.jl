#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-07: add CanopyOptics struct (will be a field for canopy structure)
#     2022-Jun-07: add more fields: fo, fs, po, ps, pso, _Co, _Cs, _So, _Ss, _abs_fo, _abs_fs, _abs_fs_fo, _cos_θ_azi_raa, _fs_fo, _tmp_mat_incl_azi_1, _tmp_mat_incl_azi_2
#     2022-Jun-08: add more fields: ρ_dd, ρ_sd, σ_ddb, σ_ddf, σ_dob, σ_dof, σ_sdb, σ_sdf, σ_so, τ_dd, τ_sd, _tmp_vec_λ, _ρ_dd, _ρ_sd, _τ_dd, _τ_sd
#     2022-Jun-09: rename variables to be more descriptive
#     2022-Jun-10: add more fields: p_sunlit, ϵ, ρ_lw, τ_lw, _mat_down, _mat_down, _tmp_vec_azi, _ρ_lw, _τ_lw
#     2022-Jun-10: add more fields for sif calculations
#     2022-Jun-13: add more fields for sif calculations
#     2022-Jun-13: remove unnecessary cache variables
#     2022-Jun-15: rename to HyperspectralMLCanopyOpticalProperty
#     2022-Jul-19: use kwdef for the constructor
#     2022-Jul-19: add dimension control to struct
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Verhoef LIDF algorithm

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralMLCanopyOpticalProperty{FT<:AbstractFloat}
    # Dimensions
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Dimension of canopy layers"
    DIM_LAYER::Int = 20
    "Dimension of SIF wave length bins"
    DIM_SIF::Int = 29
    "Dimension of SIF excitation wave length bins"
    DIM_SIFE::Int = 45
    "Dimension of short wave length bins"
    DIM_WL::Int = 114

    # Diagnostic variables
    "Backward diffuse->diffuse scatter weight"
    ddb::FT = 0
    "Forward diffuse->diffuse scatter weight"
    ddf::FT = 0
    "Backward diffuse->observer scatter weight"
    dob::FT = 0
    "Forward diffuse->observer scatter weight"
    dof::FT = 0
    "Conversion factor fo for angle towards observer at different inclination and azimuth angles"
    fo::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Conversion factor fs for angles from solar at different inclination and azimuth angles"
    fs::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Observer direction beam extinction coefficient weight (diffuse)"
    ko::FT = 0
    "Solar direction beam extinction coefficient weight (direct)"
    ks::FT = 0
    "Probability of directly viewing a leaf in solar direction at different layers"
    p_sunlit::Vector{FT} = zeros(FT, DIM_LAYER)
    "Probability of directly viewing a leaf in observer direction at different layer boundaries"
    po::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Probability of directly viewing a leaf in solar direction at different layer boundaries"
    ps::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Bi-directional probability of directly viewing a leaf at different layer boundaries (solar->canopy->observer)"
    pso::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Directional->diffuse backscatter weight"
    sdb::FT = 0
    "Directional->diffuse forward scatter weight"
    sdf::FT = 0
    "Solar directional->observer weight of specular2directional backscatter coefficient"
    sob::FT = 0
    "Solar directional->observer weight of specular2directional forward coefficient"
    sof::FT = 0
    "Effective emissivity for different layers"
    ϵ::Vector{FT} = zeros(FT, DIM_LAYER)
    "Effective reflectance for diffuse->diffuse"
    ρ_dd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Effective reflectance for longwave radiation"
    ρ_lw::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Effective reflectance for directional->diffuse"
    ρ_sd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Backward scattering coefficient for diffuse->diffuse at different layers and wavelength bins"
    σ_ddb::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Forward scattering coefficient for diffuse->diffuse at different layers and wavelength bins"
    σ_ddf::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Backward scattering coefficient for diffuse->observer at different layers and wavelength bins"
    σ_dob::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Forward scattering coefficient for diffuse->observer at different layers and wavelength bins"
    σ_dof::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Backward scattering coefficient for solar directional->diffuse at different layers and wavelength bins"
    σ_sdb::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Forward scattering coefficient for solar directional->diffuse at different layers and wavelength bins"
    σ_sdf::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Bidirectional from solar to observer scattering coefficient at different layers and wavelength bins"
    σ_so::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Effective tranmittance for diffuse->diffuse"
    τ_dd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Effective tranmittance for longwave radiation"
    τ_lw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Effective tranmittance for solar directional->diffuse"
    τ_sd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)

    # Cache variables
    "cos(inclination) * cos(vza) at different inclination angles"
    _Co::Vector{FT} = zeros(FT, DIM_INCL)
    "cos(inclination) * cos(sza) at different inclination angles"
    _Cs::Vector{FT} = zeros(FT, DIM_INCL)
    "sin(inclination) * sin(vza) at different inclination angles"
    _So::Vector{FT} = zeros(FT, DIM_INCL)
    "sin(inclination) * sin(sza) at different inclination angles"
    _Ss::Vector{FT} = zeros(FT, DIM_INCL)
    "abs of fo"
    _abs_fo::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "abs of fs"
    _abs_fs::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "abs of fs * fo"
    _abs_fs_fo::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Weighted sum of cos²(inclination)"
    _bf::FT = 0
    "Cosine of Θ_AZI - raa"
    _cos_θ_azi_raa::Vector{FT} = zeros(FT, DIM_AZI)
    "fo * cos Θ_INCL"
    _fo_cos_θ_incl::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "fs * cos Θ_INCL"
    _fs_cos_θ_incl::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "fs * fo"
    _fs_fo::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Outgoing beam extinction coefficient weights at different inclination angles"
    _ko::Vector{FT} = zeros(FT, DIM_INCL)
    "Solar beam extinction coefficient weights at different inclination angles"
    _ks::Vector{FT} = zeros(FT, DIM_INCL)
    "Upwelling matrix for SIF excitation"
    _mat⁺::Matrix{FT} = zeros(FT, DIM_SIF, DIM_SIFE)
    "Downwelling matrix for SIF excitation"
    _mat⁻::Matrix{FT} = zeros(FT, DIM_SIF, DIM_SIFE)
    "Backward scattering coefficients at different inclination angles"
    _sb::Vector{FT} = zeros(FT, DIM_INCL)
    "Forward scattering coefficients at different inclination angles"
    _sf::Vector{FT} = zeros(FT, DIM_INCL)
    "Temporary cache used for matrix adding up purpose (DIM_INCL * DIM_AZI)"
    _tmp_mat_incl_azi_1::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Temporary cache used for matrix adding up purpose (DIM_INCL * DIM_AZI)"
    _tmp_mat_incl_azi_2::Matrix{FT} = zeros(FT, DIM_INCL, DIM_AZI)
    "Temporary cache used for vector operations (DIM_AZI)"
    _tmp_vec_azi::Vector{FT} = zeros(FT, DIM_AZI)
    "Temporary cache used for vector operations (DIM_LAYER)"
    _tmp_vec_layer::Vector{FT} = zeros(FT, DIM_LAYER)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_1::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_2::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_3::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_4::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_5::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF information"
    _tmp_vec_sif_6::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_1::Vector{FT} = zeros(FT, DIM_SIFE)
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_2::Vector{FT} = zeros(FT, DIM_SIFE)
    "Cache variable to store the SIF excitation information"
    _tmp_vec_sife_3::Vector{FT} = zeros(FT, DIM_SIFE)
    "Temporary cache used for vector operations (DIM_WL)"
    _tmp_vec_λ::Vector{FT} = zeros(FT, DIM_WL)
    "Reflectance for diffuse->diffuse at each canopy layer"
    _ρ_dd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Reflectance for longwave radiation at each canopy layer"
    _ρ_lw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Reflectance for solar directional->diffuse at each canopy layer"
    _ρ_sd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Tranmittance for diffuse->diffuse at each canopy layer"
    _τ_dd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Tranmittance for longwave radiation at each canopy layer"
    _τ_lw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Tranmittance for solar directional->diffuse at each canopy layer"
    _τ_sd::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Tranmittance for solar directional->directional at each canopy layer"
    _τ_ss::Vector{FT} = zeros(FT, DIM_LAYER)
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-15: add abstract type for canopy radiation profile
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractCanopy:
- [`BroadbandSLCanopyRadiationProfile`](@ref)
- [`HyperspectralMLCanopyRadiationProfile`](@ref)

"""
abstract type AbstractCanopyRadiationProfile{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-15: add struct for broadband radiation
#     2022-Jun-16: add cache values for diffuse and direct radiation
#     2022-Jun-16: add more variable to store partitions and radiations
#     2022-Jul-19: use kwdef for the constructor
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to store canopy radiation profiles

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandSLCanopyRadiationProfile{FT<:AbstractFloat} <: AbstractCanopyRadiationProfile{FT}
    # dimensions
    "Dimension of inclination angles"
    DIM_INCL::Int = 9

    # diagnostic variables that change with time
    "Mean shaded leaf APAR (per leaf area) in μmol m⁻² s⁻¹"
    apar_shaded::FT = 0
    "Mean sunlit leaf APAR (per leaf area) in μmol m⁻² s⁻¹"
    apar_sunlit::FT = 0
    "Weighted extinction coefficient for diffuse radiation (ratio between projected area to true leaf area)"
    k_diffuse::FT = 0
    "Weighted extinction coefficient for direct radiation (ratio between projected area to true leaf area)"
    k_direct::FT = 0
    "Total shaded leaf area index"
    lai_shaded::FT = 0
    "Total sunlit leaf area index"
    lai_sunlit::FT = 0
    "Mean shaded leaf PAR (per leaf area) in μmol m⁻² s⁻¹"
    par_shaded::FT = 0
    "Mean sunlit leaf PAR (per leaf area) in μmol m⁻² s⁻¹"
    par_sunlit::FT = 0
    "Net absorbed radiation for shaded leaves `[W m⁻²]`"
    r_net_shaded::FT = 0
    "Net absorbed radiation for sunlit leaves `[W m⁻²]`"
    r_net_sunlit::FT = 0

    # caches to speed up calculations
    "Extinction coefficient for diffuse radiation at different leaf inclination angles"
    _k_diffuse::Vector{FT} = zeros(FT, DIM_INCL)
    "Extinction coefficient for direct radiation at different leaf inclination angles"
    _k_direct::Vector{FT} = zeros(FT, DIM_INCL)
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-09: migrate CanopyRads as HyperspectralMLCanopyRadiationProfile
#     2022-Jun-09: add fields: albedo, apar_shaded, apar_sunlit, e_net_diffuse, e_net_direct, e_o, e_v, par_shaded, par_sunlit, r_net
#     2022-Jun-10: add fields: e_sum_diffuse, e_sum_direct, par_in, par_in_diffuse, par_in_direct, par_shaded, par_sunlit, _par_shaded, _par_sunlit
#     2022-Jun-10: add fields: r_net_sw, r_net_sw_shaded, r_net_sw_sunlit, r_lw, r_lw_down, r_lw_up, _r_emit_down, _r_emit_up
#     2022-Jun-10: add more fields for SIF
#     2022-Jun-13: add more fields for sif calculations
#     2022-Jun-15: rename to HyperspectralMLCanopyRadiationProfile
#     2022-Jun-27: move ppar_sunlit and ppar_shaded to Leaves2D
#     2022-Jul-19: add dimension control to struct
#     2022-Aug-30: rename sif_obs_ssoil to sif_obs_soil (typo fix, non-breaking)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to store canopy radiation profiles

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralMLCanopyRadiationProfile{FT<:AbstractFloat} <: AbstractCanopyRadiationProfile{FT}
    # Dimensions
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Dimension of canopy layers"
    DIM_LAYER::Int = 20
    "Dimension of PAR wave length bins"
    DIM_PAR::Int = 35
    "Dimension of SIF wave length bins"
    DIM_SIF::Int = 29
    "Dimension of short wave length bins"
    DIM_WL::Int = 114

    # Diagnostic variables
    "Albedo towards the viewing direction"
    albedo::Vector{FT} = zeros(FT, DIM_WL)
    "Mean APAR for shaded leaves `[μmol m⁻² s⁻¹]`"
    apar_shaded::Vector{FT} = zeros(FT, DIM_LAYER)
    "APAR for sunlit leaves `[μmol m⁻² s⁻¹]`"
    apar_sunlit::Array{FT,3} = zeros(FT, DIM_INCL, DIM_AZI, DIM_LAYER)
    "Downwelling diffuse short-wave radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_diffuse_down::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Upwelling diffuse short-wave radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_diffuse_up::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Solar directly radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_direct::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Net diffuse radiation at each canopy layer for APAR `[mW m⁻² nm⁻¹]`"
    e_net_diffuse::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Net direct radiation at each canopy layer for APAR `[mW m⁻² nm⁻¹]`"
    e_net_direct::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Total radiation towards the viewing direction `[mW m⁻² nm⁻¹]`"
    e_o::Vector{FT} = zeros(FT, DIM_WL)
    "Sum diffuse radiation at each canopy layer for PAR `[mW m⁻² nm⁻¹]`"
    e_sum_diffuse::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Sum direct radiation at each canopy layer for PAR `[mW m⁻² nm⁻¹]`"
    e_sum_direct::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER)
    "Radiation towards the viewing direction per layer (including soil) `[mW m⁻² nm⁻¹]`"
    e_v::Matrix{FT} = zeros(FT, DIM_WL, DIM_LAYER+1)
    "Total incoming radiation PAR `[μmol m⁻² s⁻¹]`"
    par_in::FT = 0
    "Diffuse incoming radiation PAR `[μmol m⁻² s⁻¹]`"
    par_in_diffuse::FT = 0
    "Direct incoming radiation PAR `[μmol m⁻² s⁻¹]`"
    par_in_direct::FT = 0
    "Mean PAR for shaded leaves (before absorption) `[μmol m⁻² s⁻¹]`"
    par_shaded::Vector{FT} = zeros(FT, DIM_LAYER)
    "PAR for sunlit leaves (before absorption) `[μmol m⁻² s⁻¹]`"
    par_sunlit::Array{FT,3} = zeros(FT, DIM_INCL, DIM_AZI, DIM_LAYER)
    "Longwave energy flux from leaves per leaf area (one side) `[W m⁻²]`"
    r_lw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Downwelling longwave energy flux `[W m⁻²]`"
    r_lw_down::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Upwelling longwave energy flux `[W m⁻²]`"
    r_lw_up::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Net longwave energy absorption for all leaves `[W m⁻²]`"
    r_net_lw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Net shortwave energy absorption for all leaves `[W m⁻²]`"
    r_net_sw::Vector{FT} = zeros(FT, DIM_LAYER)
    "Net shortwave energy absorption for shaded leaves `[W m⁻²]`"
    r_net_sw_shaded::Vector{FT} = zeros(FT, DIM_LAYER)
    "Net shortwave energy absorption for sunlit leaves `[W m⁻²]`"
    r_net_sw_sunlit::Vector{FT} = zeros(FT, DIM_LAYER)
    "Downwelling SIF for sunlit leaves at each wavelength for a layer"
    s_layer_down::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
    "Upwelling SIF for sunlit leaves at each wavelength for a layer"
    s_layer_up::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
    "Downwelling SIF"
    sif_down::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER+1)
    "SIF at observer direction"
    sif_obs::Vector{FT} = zeros(FT, DIM_SIF)
    "SIF at observer direction from shaded APAR"
    sif_obs_shaded::Vector{FT} = zeros(FT, DIM_SIF)
    "SIF at observer direction from scattering"
    sif_obs_scatter::Vector{FT} = zeros(FT, DIM_SIF)
    "SIF at observer direction from soil reflection"
    sif_obs_soil::Vector{FT} = zeros(FT, DIM_SIF)
    "SIF at observer direction from sunlit APAR"
    sif_obs_sunlit::Vector{FT} = zeros(FT, DIM_SIF)
    "Upwelling SIF"
    sif_up::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER+1)

    # Cache variables
    "Mean APAR for shaded leaves per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _apar_shaded::Vector{FT} = zeros(FT, DIM_PAR)
    "APAR for sunlit leaves per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _apar_sunlit::Vector{FT} = zeros(FT, DIM_PAR)
    "Mean PAR for shaded leaves per wavelength (before absorption) `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _par_shaded::Vector{FT} = zeros(FT, DIM_PAR)
    "PAR for sunlit leaves per wavelength (before absorption) `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _par_sunlit::Vector{FT} = zeros(FT, DIM_PAR)
    "Mean APAR for shaded leaves for photosynthesis per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _ppar_shaded::Vector{FT} = zeros(FT, DIM_PAR)
    "APAR for sunlit leaves for photosynthesis per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _ppar_sunlit::Vector{FT} = zeros(FT, DIM_PAR)
    "Downwelling longwave energy flux cache `[W m⁻²]`"
    _r_emit_down::Vector{FT} = zeros(FT, DIM_LAYER)
    "Upwelling longwave energy flux cache `[W m⁻²]`"
    _r_emit_up::Vector{FT} = zeros(FT, DIM_LAYER+1)
    "Downwelling SIF for sunlit leaves at each wavelength"
    _s_emit_down::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
    "Upwelling SIF for sunlit leaves at each wavelength"
    _s_emit_up::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER+1)
    "Downwelling SIF for shaded leaves at each wavelength"
    _s_shaded_down::Vector{FT} = zeros(FT, DIM_SIF)
    "Upwelling SIF for shaded leaves at each wavelength"
    _s_shaded_up::Vector{FT} = zeros(FT, DIM_SIF)
    "Downwelling SIF for sunlit leaves at each wavelength"
    _s_sunlit_down::Vector{FT} = zeros(FT, DIM_SIF)
    "Upwelling SIF for sunlit leaves at each wavelength"
    _s_sunlit_up::Vector{FT} = zeros(FT, DIM_SIF)
    "Cache to compute SIF at observer direction from shaded APAR"
    _sif_obs_shaded::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
    "Cache to compute SIF at observer direction from scattering"
    _sif_obs_scatter::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
    "Cache to compute SIF at observer direction from sunlit APAR"
    _sif_obs_sunlit::Matrix{FT} = zeros(FT, DIM_SIF, DIM_LAYER)
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-02: add abstract type for LIDF algorithms
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractLIDFAlgorithm:
- [`BetaLIDF`](@ref)
- [`VerhoefLIDF`](@ref)

"""
abstract type AbstractLIDFAlgorithm{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2023-May-22: add beta function method
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Beta LIDF algorithm

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BetaLIDF{FT<:AbstractFloat} <: AbstractLIDFAlgorithm{FT}
    # General model information
    "Leaf inclination angle distribution function parameter a"
    A::FT = 1
    "Leaf inclination angle distribution function parameter b"
    B::FT = 1
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-02: migrate from CanopyLayers
#     2022-Jun-02: rename Canopy4RT to HyperspectralMLCanopy
#     2022-Jun-02: abstractize LIDF as a field
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Verhoef LIDF algorithm

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct VerhoefLIDF{FT<:AbstractFloat} <: AbstractLIDFAlgorithm{FT}
    # General model information
    "Leaf inclination angle distribution function parameter a"
    A::FT = 0
    "Leaf inclination angle distribution function parameter b"
    B::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-02: add abstract type for canopy structure
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractCanopy:
- [`BroadbandSLCanopy`](@ref)
- [`HyperspectralMLCanopy`](@ref)

"""
abstract type AbstractCanopy{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-15: add struct for broadband radiative transfer scheme such as two leaf model
#     2022-Jun-15: add more cache variables
#     2022-Jun-15: add radiation profile
#     2022-Jun-15: remove RATIO_HV to compute the coefficient numerically
#     2022-Jun-16: remove some cache variables
#     2022-Jun-16: add fields: Θ_INCL_BNDS
#     2023-May-22: add sypport to BetaLIDF
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save single layer broadband canopy parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandSLCanopy{FT<:AbstractFloat} <: AbstractCanopy{FT}
    # Dimensions
    "Dimension of inclination angles"
    DIM_INCL::Int = 9

    # Embedded structures
    "Leaf inclination angle distribution function algorithm"
    LIDF::Union{BetaLIDF{FT}, VerhoefLIDF{FT}} = VerhoefLIDF{FT}()
    "Canopy radiation profiles"
    RADIATION::BroadbandSLCanopyRadiationProfile{FT} = BroadbandSLCanopyRadiationProfile{FT}(DIM_INCL = DIM_INCL)

    # Geometry information
    "Inclination angle distribution"
    P_INCL::Vector{FT} = ones(FT, DIM_INCL) ./ DIM_INCL
    "Bounds of inclination angles `[°]`"
    Θ_INCL_BNDS::Matrix{FT} = FT[ collect(FT, range(0, 90; length=DIM_INCL+1))[1:end-1] collect(FT, range(0, 90; length=DIM_INCL+1))[2:end] ]
    "Mean inclination angles `[°]`"
    Θ_INCL::Vector{FT} = [ (Θ_INCL_BNDS[_i,1] + Θ_INCL_BNDS[_i,2]) / 2 for _i in 1:DIM_INCL ]

    # Prognostic variables
    "Clumping index"
    ci::FT = 1
    "Leaf area index"
    lai::FT = 3
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-02: migrate from CanopyLayers
#     2022-Jun-02: rename Canopy4RT to HyperspectralMLCanopy
#     2022-Jun-02: abstractize LIDF as a field
#     2022-Jun-07: add cache variable _1_AZI, _COS²_Θ_INCL, _COS_Θ_INCL_AZI, _COS²_Θ_INCL_AZI
#     2022-Jun-07: remove cache variable _cos_θ_azi_raa, _vol_scatter
#     2022-Jun-09: add new field: APAR_CAR, RADIATION, WLSET
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jun-15: rename to HyperspectralMLCanopyOpticalProperty and HyperspectralMLCanopyRadiationProfile
#     2022-Jun-16: remove some cache variables
#     2022-Jul-22: remove field APAR_CAR
#     2022-Aug-30: add field LHA (moved from spac)
#     2023-May-22: add sypport to BetaLIDF
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save multiple layer hyperspectral canopy parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralMLCanopy{FT<:AbstractFloat} <: AbstractCanopy{FT}
    # Dimensions
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Dimension of canopy layers"
    DIM_LAYER::Int = 20

    # General model information
    "Hot spot parameter"
    HOT_SPOT::FT = 0.05

    # Embedded structures
    "Hyperspectral absorption features of different leaf components"
    LHA::HyperspectralAbsorption{FT} = HyperspectralAbsorption{FT}()
    "Leaf inclination angle distribution function algorithm"
    LIDF::Union{BetaLIDF{FT}, VerhoefLIDF{FT}} = VerhoefLIDF{FT}()
    "Wave length set used to paramertize other variables"
    WLSET::WaveLengthSet{FT} = WaveLengthSet{FT}()
    "Canopy optical properties"
    OPTICS::HyperspectralMLCanopyOpticalProperty{FT} = HyperspectralMLCanopyOpticalProperty{FT}(
                DIM_AZI = DIM_AZI,
                DIM_INCL = DIM_INCL,
                DIM_LAYER = DIM_LAYER,
                DIM_SIF = WLSET.DIM_SIF,
                DIM_SIFE = WLSET.DIM_SIFE,
                DIM_WL = WLSET.DIM_WL)
    "Canopy radiation profiles"
    RADIATION::HyperspectralMLCanopyRadiationProfile{FT} = HyperspectralMLCanopyRadiationProfile{FT}(
                DIM_AZI = DIM_AZI,
                DIM_INCL = DIM_INCL,
                DIM_LAYER = DIM_LAYER,
                DIM_PAR = WLSET.DIM_PAR,
                DIM_SIF = WLSET.DIM_SIF,
                DIM_WL = WLSET.DIM_WL)

    # Geometry information
    "Inclination angle distribution"
    P_INCL::Vector{FT} = ones(FT, DIM_INCL) ./ DIM_INCL
    "Mean azimuth angles `[°]`"
    Θ_AZI::Vector{FT} = collect(FT, range(0, 360; length=DIM_AZI+1))[1:end-1] .+ 360 / DIM_AZI / 2
    "Bounds of inclination angles `[°]`"
    Θ_INCL_BNDS::Matrix{FT} = FT[ collect(FT, range(0, 90; length=DIM_INCL+1))[1:end-1] collect(FT, range(0, 90; length=DIM_INCL+1))[2:end] ]
    "Mean inclination angles `[°]`"
    Θ_INCL::Vector{FT} = [ (Θ_INCL_BNDS[_i,1] + Θ_INCL_BNDS[_i,2]) / 2 for _i in 1:DIM_INCL ]
    "Clumping structure a"
    Ω_A::FT = 1
    "Clumping structure b"
    Ω_B::FT = 0

    # Prognostic variables
    "Clumping index"
    ci::FT = 1
    "Leaf area index"
    lai::FT = 3
    "Leaf area index distribution"
    δlai::Vector{FT} = lai .* ones(FT, DIM_LAYER) ./ DIM_LAYER

    # Cache variables
    "Ones with the length of Θ_AZI"
    _1_AZI::Vector{FT} = ones(FT, DIM_AZI)
    "Cosine of Θ_AZI"
    _COS_Θ_AZI::Vector{FT} = cosd.(Θ_AZI)
    "Square of cosine of Θ_INCL"
    _COS²_Θ_INCL::Vector{FT} = cosd.(Θ_INCL) .^ 2
    "Square of cosine of Θ_INCL at different azimuth angles"
    _COS²_Θ_INCL_AZI::Matrix{FT} = (cosd.(Θ_INCL) .^ 2) * _1_AZI'
    "Cache for level boundary locations"
    _x_bnds::Vector{FT} = [0; [sum(δlai[1:_i]) for _i in 1:DIM_LAYER]] ./ -lai
end
