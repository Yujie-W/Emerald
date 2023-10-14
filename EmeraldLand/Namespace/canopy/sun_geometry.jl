# This file contains structs for sun geometry (apart from sensor geometry)

#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-09: add struct SunGeometryState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the state variables of the sun geometry.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SunGeometryState{FT}
    "Solar azimuth angle `[°]`, a function of time"
    saa::FT = 180
    "Solar zenith angle `[°]`, a function of lat and time"
    sza::FT = 30
end;


#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-09: add struct SunGeometryAuxil
#     2023-Oct-13: add field albedo
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the auxiliary variables of the sun geometry.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SunGeometryAuxil{FT}
    # Scattering coefficients
    "Backward direct->diffuse scatter weight"
    sdb::FT = 0
    "Forward direct->diffuse scatter weight"
    sdf::FT = 0

    # Extinction coefficient related
    "Solar direction beam extinction coefficient weight (direct)"
    ks::FT = 0
    "Probability of directly viewing a leaf in solar direction at different layers"
    p_sunlit::Vector{FT}

    # Extinction coefficient related (for different inclination angles)
    "cos(inclination) * cos(sza) at different inclination angles"
    Cs_incl::Vector{FT}
    "sin(inclination) * sin(sza) at different inclination angles"
    Ss_incl::Vector{FT}
    "Solar beam extinction coefficient weights at different inclination angles"
    ks_incl::Vector{FT}
    "Cs >= Ss ? FT(π) : acos(-Cs/Ss)"
    βs_incl::Vector{FT}

    # Matrix used for solar radiation
    "Conversion factor fs for angles from solar at different inclination and azimuth angles"
    fs::Matrix{FT}
    "Absolute value of fs"
    fs_abs::Matrix{FT}
    "Mean fs_abs at different azimuth angles"
    fs_abs_mean::Vector{FT}
    "fs * cos Θ_INCL"
    fs_cos²_incl::Matrix{FT}

    # Scattering coefficients per leaf area
    "Backward scattering coefficient for solar directional->diffuse at different layers and wavelength bins"
    sdb_leaf::Matrix{FT}
    "Forward scattering coefficient for solar directional->diffuse at different layers and wavelength bins"
    sdf_leaf::Matrix{FT}

    # Reflectance and tranmittance per canopy layer (no denominator correction made yet)
    "Reflectance for solar directional->diffuse at each canopy layer"
    ρ_sd_layer::Matrix{FT}
    "Tranmittance for solar directional->diffuse at each canopy layer"
    τ_sd_layer::Matrix{FT}
    "Tranmittance for solar directional->directional at each canopy layer (wavelength independent)"
    τ_ss_layer::Vector{FT}

    # Effective reflectance and tranmittance per canopy layer (including the denominator correction)
    "Effective reflectance for directional->diffuse"
    ρ_sd::Matrix{FT}
    "Effective tranmittance for solar directional->diffuse"
    τ_sd::Matrix{FT}

    # Canopy radiation profiles
    "Albedo of the soil-canopy of the entire hemisphere"
    albedo::Vector{FT}
    "Downwelling diffuse short-wave radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_difꜜ::Matrix{FT}
    "Upwelling diffuse short-wave radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_difꜛ::Matrix{FT}
    "Solar directly radiation at each canopy layer boundary `[mW m⁻² nm⁻¹]`"
    e_dirꜜ::Matrix{FT}

    # Radiation at each canopy layer used for PAR, APAR, PPAR calculations
    "Net diffuse radiation at each canopy layer for APAR `[mW m⁻² nm⁻¹]`"
    e_net_dif::Matrix{FT}
    "Net direct radiation at each canopy layer for APAR `[mW m⁻² nm⁻¹]`"
    e_net_dir::Matrix{FT}

    # Net radiation used for energy budget
    "Net shortwave energy absorption for all leaves per leaf area `[W m⁻²]`"
    r_net_sw::Vector{FT}

    # Fluorescence related fluxes (not computing the fluorescence spectrum at sensor direction)
    "SIF emission per layer at chloroplast level"
    e_sif_chl::Matrix{FT}
    "Downward SIF emitted per layer"
    e_sifꜜ_layer::Matrix{FT}
    "Upward SIF emitted per layer"
    e_sifꜛ_layer::Matrix{FT}
    "Downward effective emitted SIF per layer"
    e_sifꜜ_emit::Matrix{FT}
    "Upward effective emitted SIF per layer"
    e_sifꜛ_emit::Matrix{FT}
    "Downward SIF per layer"
    e_sifꜜ::Matrix{FT}
    "Upward SIF per layer"
    e_sifꜛ::Matrix{FT}

    # PAR, APAR, PPAR flux density in each layer
    "Mean APAR for shaded leaves per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _apar_shaded::Vector{FT}
    "APAR for sunlit leaves per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _apar_sunlit::Vector{FT}
    "Mean APAR for shaded leaves for photosynthesis per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _ppar_shaded::Vector{FT}
    "APAR for sunlit leaves for photosynthesis per wavelength `[μmol m⁻² s⁻¹ nm⁻¹]`"
    _ppar_sunlit::Vector{FT}

    # Fluorescence related cache variables (related to chlorophyll level fluorescence)
    "Diffuse radiation that can excite SIF emission per wavelength"
    _e_dif_sife::Vector{FT}
    "Direct radiation that can excite SIF emission per wavelength"
    _e_dir_sife::Vector{FT}
    "SIF emission from diffuse radiation per wavelength"
    _e_dif_sif::Vector{FT}
    "SIF emission from direct radiation per wavelength"
    _e_dir_sif::Vector{FT}
    "SIF emission from shaded leaves per wavelength from diffuse radiation"
    _e_dif_shaded::Vector{FT}
    "SIF emission from sunlit leaves per wavelength from diffuse radiation"
    _e_dif_sunlit::Vector{FT}
    "SIF emission from sunlit leaves per wavelength from direct radiation"
    _e_dir_sunlit::Vector{FT}

    # Fluorescence related cache variables (related to leaf level fluorescence after leaf reabsorption)
    "Downward diffuse radiation that can excite SIF emission per wavelength"
    _e_difꜜ_sife::Vector{FT}
    "Upward diffuse radiation that can excite SIF emission per wavelength"
    _e_difꜛ_sife::Vector{FT}
    "Downward direct radiation that can excite SIF emission per wavelength"
    _e_dirꜜ_sife::Vector{FT}
    "SIF emission from downward diffuse radiation per wavelength (mean)"
    _e_difꜜ_sif_mean::Vector{FT}
    "SIF emission from downward diffuse radiation per wavelength (diff)"
    _e_difꜜ_sif_diff::Vector{FT}
    "SIF emission from upward diffuse radiation per wavelength (mean)"
    _e_difꜛ_sif_mean::Vector{FT}
    "SIF emission from upward diffuse radiation per wavelength (diff)"
    _e_difꜛ_sif_diff::Vector{FT}
    "SIF emission from downward direct radiation per wavelength (mean)"
    _e_dirꜜ_sif_mean::Vector{FT}
    "SIF emission from downward direct radiation per wavelength (diff)"
    _e_dirꜜ_sif_diff::Vector{FT}
    "Downward SIF emissions from shaded leaves"
    _sif_shadedꜜ::Vector{FT}
    "Upward SIF emissions from shaded leaves"
    _sif_shadedꜛ::Vector{FT}
    "Downward SIF emissions from sunlit leaves"
    _sif_sunlitꜜ::Vector{FT}
    "Upward SIF emissions from sunlit leaves"
    _sif_sunlitꜛ::Vector{FT}

    # cache variables
    "A temporary matrix with size of DIM_INCL and DIM_AZI"
    _mat_incl_azi::Matrix{FT}
    "A temporary vector with length of DIM_AZI"
    _vec_azi::Vector{FT}
end;

SunGeometryAuxil(config::SPACConfiguration{FT}) where {FT} = SunGeometryAuxil{FT}(
            p_sunlit         = zeros(FT, config.DIM_LAYER),
            Cs_incl          = zeros(FT, config.DIM_INCL),
            Ss_incl          = zeros(FT, config.DIM_INCL),
            ks_incl          = zeros(FT, config.DIM_INCL),
            βs_incl          = zeros(FT, config.DIM_INCL),
            fs               = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fs_abs           = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            fs_abs_mean      = zeros(FT, config.DIM_AZI),
            fs_cos²_incl     = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            sdb_leaf         = zeros(FT, config.DIM_WL, config.DIM_LAYER),
            sdf_leaf         = zeros(FT, config.DIM_WL, config.DIM_LAYER),
            ρ_sd_layer       = zeros(FT, config.DIM_WL, config.DIM_LAYER),
            τ_sd_layer       = zeros(FT, config.DIM_WL, config.DIM_LAYER),
            τ_ss_layer       = zeros(FT, config.DIM_LAYER),
            ρ_sd             = zeros(FT, config.DIM_WL, config.DIM_LAYER + 1),
            τ_sd             = zeros(FT, config.DIM_WL, config.DIM_LAYER),
            albedo           = zeros(FT, config.DIM_WL),
            e_difꜜ           = zeros(FT, config.DIM_WL, config.DIM_LAYER + 1),
            e_difꜛ           = zeros(FT, config.DIM_WL, config.DIM_LAYER + 1),
            e_dirꜜ           = zeros(FT, config.DIM_WL, config.DIM_LAYER + 1),
            e_net_dif        = zeros(FT, config.DIM_WL, config.DIM_LAYER),
            e_net_dir        = zeros(FT, config.DIM_WL, config.DIM_LAYER),
            r_net_sw         = zeros(FT, config.DIM_LAYER),
            e_sif_chl        = zeros(FT, config.DIM_SIF, config.DIM_LAYER),
            e_sifꜜ_layer     = zeros(FT, config.DIM_SIF, config.DIM_LAYER),
            e_sifꜛ_layer     = zeros(FT, config.DIM_SIF, config.DIM_LAYER),
            e_sifꜜ_emit      = zeros(FT, config.DIM_SIF, config.DIM_LAYER),
            e_sifꜛ_emit      = zeros(FT, config.DIM_SIF, config.DIM_LAYER + 1),
            e_sifꜜ           = zeros(FT, config.DIM_SIF, config.DIM_LAYER + 1),
            e_sifꜛ           = zeros(FT, config.DIM_SIF, config.DIM_LAYER + 1),
            _apar_shaded     = zeros(FT, config.DIM_PAR),
            _apar_sunlit     = zeros(FT, config.DIM_PAR),
            _ppar_shaded     = zeros(FT, config.DIM_PAR),
            _ppar_sunlit     = zeros(FT, config.DIM_PAR),
            _e_dif_sife      = zeros(FT, config.DIM_SIFE),
            _e_dir_sife      = zeros(FT, config.DIM_SIFE),
            _e_dif_sif       = zeros(FT, config.DIM_SIF),
            _e_dir_sif       = zeros(FT, config.DIM_SIF),
            _e_dif_shaded    = zeros(FT, config.DIM_SIF),
            _e_dif_sunlit    = zeros(FT, config.DIM_SIF),
            _e_dir_sunlit    = zeros(FT, config.DIM_SIF),
            _e_difꜜ_sife     = zeros(FT, config.DIM_SIFE),
            _e_difꜛ_sife     = zeros(FT, config.DIM_SIFE),
            _e_dirꜜ_sife     = zeros(FT, config.DIM_SIFE),
            _e_difꜜ_sif_mean = zeros(FT, config.DIM_SIF),
            _e_difꜜ_sif_diff = zeros(FT, config.DIM_SIF),
            _e_difꜛ_sif_mean = zeros(FT, config.DIM_SIF),
            _e_difꜛ_sif_diff = zeros(FT, config.DIM_SIF),
            _e_dirꜜ_sif_mean = zeros(FT, config.DIM_SIF),
            _e_dirꜜ_sif_diff = zeros(FT, config.DIM_SIF),
            _sif_shadedꜜ     = zeros(FT, config.DIM_SIF),
            _sif_shadedꜛ     = zeros(FT, config.DIM_SIF),
            _sif_sunlitꜜ     = zeros(FT, config.DIM_SIF),
            _sif_sunlitꜛ     = zeros(FT, config.DIM_SIF),
            _mat_incl_azi    = zeros(FT, config.DIM_INCL, config.DIM_AZI),
            _vec_azi         = zeros(FT, config.DIM_AZI),
);


#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-09: add struct SunGeometry
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to store the sun geometry.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SunGeometry{FT}
    "State variables"
    state::SunGeometryState{FT} = SunGeometryState{FT}()
    "Auxiliary variables"
    auxil::SunGeometryAuxil{FT}
end;

SunGeometry(config::SPACConfiguration{FT}) where {FT} = SunGeometry{FT}(auxil = SunGeometryAuxil(config));
