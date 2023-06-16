#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2021-Sep-30: add abstract type for soil vulnerability curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractSoilVC:
- [`BrooksCorey`](@ref)
- [`VanGenuchten`](@ref)

"""
abstract type AbstractSoilVC{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Sep-30: define this structure with no default constructor
#     2022-Oct-19: make struct mutable
#     2023-May-12: add BrooksCorey constructor
# Sources
#     https://ral.ucar.edu/sites/default/files/public/product-tool/unified-noah-lsm/parameters/SOILPARM.TBL
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Brooks Corey soil parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BrooksCorey{FT<:AbstractFloat} <:AbstractSoilVC{FT}
    # General model information
    "Soil b"
    B::FT
    "Maximum soil hydraulic conductivity at 25 °C `[mol m⁻¹ s⁻¹ MPa⁻¹]`"
    K_MAX::FT
    "Soil type"
    TYPE::String
    "Potential at saturation `[MPa]`"
    Ψ_SAT::FT
    "Saturated soil volumetric water content"
    Θ_SAT::FT
    "Residual soil volumetric water content"
    Θ_RES::FT
end

"""

    BrooksCorey{FT}(catg::Int) where {FT<:AbstractFloat}

Return a BrooksCorey soil VC, given
- `catg` Soil texture catergory (must be within [1,19])

"""
BrooksCorey{FT}(catg::Int) where {FT<:AbstractFloat} = (
    @assert 1 <= catg <= 19 "Soil texture catergory must be within 1 to 19!";

    return BrooksCorey{FT}(
                B     = SOIL_TEXT.BB[catg],
                K_MAX = SOIL_TEXT.SATDK[catg] / GRAVITY(FT) * 1e6 / M_H₂O(FT),
                TYPE  = SOIL_TEXT.NAME[catg],
                Ψ_SAT = SOIL_TEXT.SATPSI[catg] * ρ_H₂O(FT) * GRAVITY(FT) * 1e-6,
                Θ_SAT = SOIL_TEXT.MAXSMC[catg],
                Θ_RES = SOIL_TEXT.REFSMC[catg]
    )
);


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Sep-30: define this structure with two default constructors from an incomplete parameter set
#     2021-Sep-30: add constructor function
#     2022-Jul-13: remove a constructor method
#     2022-Oct-19: make struct mutable
#     2023-Apr-08: add another dataset for van Genuchten parameters
# Sources
#     Cosby et al. (1984) A statistical exploration of the relationships of soil moisture characteristics to the physical properties of soils
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

van Genuchten soil parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct VanGenuchten{FT<:AbstractFloat} <:AbstractSoilVC{FT}
    # General model information
    "Maximum soil hydraulic conductivity at 25 °C `[mol m⁻¹ s⁻¹ MPa⁻¹]`"
    K_MAX::FT
    "Soil n is Measure of the pore-size distribution"
    N::FT
    "Soil type"
    TYPE::String
    "Soil α is related to the inverse of the air entry suction, α > 0"
    α::FT
    "Residual soil volumetric water content"
    Θ_RES::FT
    "Saturated soil volumetric water content"
    Θ_SAT::FT

    # Parameters based on the ones above
    "Soil m = 1 - 1/n"
    M::FT = 1 - 1 / N
end


"""

    VanGenuchten{FT}(name::String) where {FT<:AbstractFloat}

Constructor for [`VanGenuchten`](@ref), given
- `name` Soil type name, need to be Sand, Loamy Sand, Sandy Loam, Loam (default), Sandy Clay Loam, Silt Loam, Silt, Clay Loam, Silty Clay Loam, Sandy Clay, Silty Clay, and Clay.

"""
VanGenuchten{FT}(name::String) where {FT<:AbstractFloat} = (
    #=
    # Parameters from Loam soil
    _p = [ 367.3476, 1.56, 0.43, 0.078, exp(-0.32) * 0.0254 / 3600];

    # switch name
    if name=="Sand"
        _p = [1479.5945, 2.68, 0.43, 0.045, exp( 0.82) * 0.0254 / 3600];
    elseif name=="Loamy Sand"
        _p = [1265.3084, 2.28, 0.41, 0.057, exp( 0.30) * 0.0254 / 3600];
    elseif name=="Sandy Loam"
        _p = [ 765.3075, 1.89, 0.41, 0.065, exp(-0.13) * 0.0254 / 3600];
    elseif name=="Loam"
        _p = [ 367.3476, 1.56, 0.43, 0.078, exp(-0.32) * 0.0254 / 3600];
    elseif name=="Sandy Clay Loam"
        _p = [ 602.0419, 1.48, 0.39, 0.100, exp(-0.20) * 0.0254 / 3600];
    elseif name=="Silt Loam"
        _p = [ 204.0820, 1.41, 0.45, 0.067, exp(-0.40) * 0.0254 / 3600];
    elseif name=="Silt"
        _p = [ 163.2656, 1.37, 0.46, 0.034, exp(-0.63) * 0.0254 / 3600];   # this k is guessed, must be within [-0.72, -0.54]
    elseif name=="Clay Loam"
        _p = [ 193.8779, 1.31, 0.41, 0.095, exp(-0.40) * 0.0254 / 3600];
    elseif name=="Silty Clay Loam"
        _p = [ 102.0410, 1.23, 0.43, 0.089, exp(-0.54) * 0.0254 / 3600];
    elseif name== "Sandy Clay"
        _p = [ 275.5107, 1.23, 0.38, 0.100, exp( 0.01) * 0.0254 / 3600];
    elseif name=="Silty Clay"
        _p = [  51.0205, 1.09, 0.36, 0.070, exp(-0.72) * 0.0254 / 3600];
    elseif name=="Clay"
        _p = [  81.6328, 1.09, 0.38, 0.068, exp(-0.86) * 0.0254 / 3600];    # K from Light clay in Cosby et al. (1984)
    else
        @warn "Soil type $(name) not recognized, use Loam instead.";
        name = "Loam";
    end;
    =#
    # https://structx.com/Soil_Properties_007.html
    # Parameters from Loam soil
    _p = [ 367.3476, 1.56, 0.43, 0.078, 7.19e-6];

    # switch name
    if name=="Sand"
        _p = [1479.5945, 2.68, 0.43, 0.045, 1.76e-4];
    elseif name=="Loamy Sand"
        _p = [1265.3084, 2.28, 0.41, 0.057, 1.56e-4];
    elseif name=="Sandy Loam"
        _p = [ 765.3075, 1.89, 0.41, 0.065, 3.45e-5];
    elseif name=="Loam"
        _p = [ 367.3476, 1.56, 0.43, 0.078, 6.94e-6];
    elseif name=="Sandy Clay Loam"
        _p = [ 602.0419, 1.48, 0.39, 0.100, 6.31e-6];
    elseif name=="Silt Loam"
        _p = [ 204.0820, 1.41, 0.45, 0.067, 7.19e-6];
    elseif name=="Silt"
        _p = [ 163.2656, 1.37, 0.46, 0.034, 7.19e-6];
    elseif name=="Clay Loam"
        _p = [ 193.8779, 1.31, 0.41, 0.095, 2.45e-6];
    elseif name=="Silty Clay Loam"
        _p = [ 102.0410, 1.23, 0.43, 0.089, 1.70e-6];
    elseif name== "Sandy Clay"
        _p = [ 275.5107, 1.23, 0.38, 0.100, 2.17e-6];
    elseif name=="Silty Clay"
        _p = [  51.0205, 1.09, 0.36, 0.070, 1.02e-6];
    elseif name=="Clay"
        _p = [  81.6328, 1.09, 0.38, 0.068, 1.28e-6];
    else
        @warn "Soil type $(name) not recognized, use Loam instead.";
        name = "Loam";
    end;

    # return a new struct
    return VanGenuchten{FT}(K_MAX = _p[5] / GRAVITY(FT) * 1e6 / M_H₂O(FT), N = _p[2], TYPE = name, α = _p[1], Θ_RES = _p[4], Θ_SAT = _p[3])
);


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jun-14: add abstract type for soil albedo
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractSoilAlbedo:
- [`BroadbandSoilAlbedo`](@ref)
- [`HyperspectralSoilAlbedo`](@ref)

"""
abstract type AbstractSoilAlbedo{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-14: add struct for broadband soil albedo
#     2022-Jun-14: make soil albedo a two-element vector for PAR and NIR
#     2022-Jul-27: add field α_CLM
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for broadband soil albedo

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BroadbandSoilAlbedo{FT<:AbstractFloat} <: AbstractSoilAlbedo{FT}
    # Constants
    "Reflectance for longwave radiation"
    ρ_LW::FT = 0.06

    # General model information
    "Whether to use CLM soil albedo scheme"
    α_CLM::Bool = false

    # Diagnostic variables
    "Net diffuse radiation at top soil `[W m⁻²]`"
    e_net_diffuse::FT = 0
    "Net direct radiation at top soil `[W m⁻²]`"
    e_net_direct::FT = 0
    "Net longwave energy absorption `[W m⁻²]`"
    r_net_lw::FT = 0
    "Net shortwave energy absorption `[W m⁻²]`"
    r_net_sw::FT = 0
    "Reflectance for shortwave radiation (for PAR and NIR)"
    ρ_sw::Vector{FT} = FT[0, 0]
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-14: add struct for hyperspectral soil albedo
#     2022-Jun-14: add constructor
#     2022-Jun-14: add fields to compute soil hyperspectral albedo in CanopyRadiativeTransfer.jl
#     2022-Jun-14: add wls in constructor function and remove n_λ
#     2022-Jul-22: rename Ρ (greek) to ρ
#     2022-Jul-27: add field α_CLM, _θ
#     2023-Jun-16: remove fields DIM_*
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for hyperspectral soil albedo

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct HyperspectralSoilAlbedo{FT<:AbstractFloat} <: AbstractSoilAlbedo{FT}
    # File path to the Netcdf dataset
    "File path to the Netcdf dataset"
    DATASET::String = LAND_2021

    # General model information
    "Whether to use CLM soil albedo scheme"
    α_CLM::Bool = false
    "Whether to fit the data from broadband to hyperspectral"
    FITTING::Bool = true

    # Constants
    "A matrix of characteristic curves"
    MAT_ρ::Matrix{FT} = FT[read_nc(DATASET, "GSV_1") read_nc(DATASET, "GSV_2") read_nc(DATASET, "GSV_3") read_nc(DATASET, "GSV_4")]
    "Reflectance for longwave radiation"
    ρ_LW::FT = 0.06

    # Diagnostic variables
    "Net diffuse radiation at top soil `[mW m⁻² nm⁻¹]`"
    e_net_diffuse::Vector{FT}
    "Net direct radiation at top soil `[mW m⁻² nm⁻¹]`"
    e_net_direct::Vector{FT}
    "Net longwave energy absorption `[W m⁻²]`"
    r_net_lw::FT = 0
    "Net shortwave energy absorption `[W m⁻²]`"
    r_net_sw::FT = 0
    "Reflectance for shortwave radiation"
    ρ_sw::Vector{FT}

    # Cache variables
    "Cache variable with length of NIR"
    _tmp_vec_nir::Vector{FT}
    "Weights of the four characteristic curves"
    _weight::Vector{FT} = zeros(FT, 4)
    "Cache variable to store ρ_PAR and ρ_NIR (a segmented curve)"
    _ρ_sw::Vector{FT}
    "Last soil moisture used to compute albedo"
    _θ::FT = -1
end

HyperspectralSoilAlbedo(config::SPACConfiguration{FT}) where {FT} = (
    (; DIM_NIR, DIM_WL) = config;

    return HyperspectralSoilAlbedo{FT}(
                e_net_diffuse = zeros(FT, DIM_WL),
                e_net_direct  = zeros(FT, DIM_WL),
                ρ_sw          = zeros(FT, DIM_WL),
                _tmp_vec_nir  = zeros(FT, DIM_NIR),
                _ρ_sw         = zeros(FT, DIM_WL),
    )
);


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2023-Jun-12: add struct SoilTraceGasses
#     2023-Jun-13: add fields for water vapor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for soil trace gasses

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilTraceGasses{FT<:AbstractFloat}
    "CH₄ mole `[mol]`"
    n_CH₄::FT = 0
    "CO₂ mole `[mol]`"
    n_CO₂::FT = 0
    "H₂O mole `[mol]`"
    n_H₂O::FT = 0
    "Nitrogen mole `[mol]`"
    n_N₂::FT = 0
    "Oxygen mole `[mol]`"
    n_O₂::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jul-13: add SoilLayer structure
#     2022-Jul-13: add field K_MAX, K_REF, k, ψ, and ∂θ∂t
#     2022-Jul-14: remove field K_REF
#     2022-Jul-14: add field ∂G∂t (renamed to ∂e∂t), ΔZ
#     2022-Jul-26: move field K_MAX to structure VanGenuchten and BrooksCorey
#     2023-Jun-13: add field TRACES for SoilTraceGasses, and ∂n∂t for the gas diffusion
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for soil layer

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilLayer{FT<:AbstractFloat}
    # Constants
    "Specific heat capacity of soil `[J K⁻¹ kg⁻¹]`"
    CP::FT = 760
    "Soil thermal conductivity `[W m⁻¹ K⁻¹]`"
    Λ_THERMAL::FT = 3
    "Dry soil density `[kg m⁻³]`"
    ρ::FT = 2650

    # Embedded structures
    "Soil trace gasses"
    TRACES::SoilTraceGasses{FT} = SoilTraceGasses{FT}()
    "Soil moisture retention curve"
    VC::Union{BrooksCorey{FT}, VanGenuchten{FT}} = VanGenuchten{FT}("Loam")

    # Geometry information
    "Depth boundaries `[m]`"
    ZS::Vector{FT} = FT[0,-1]
    "Mean depth `[m]`"
    Z::FT = (ZS[1] + ZS[2]) / 2
    "Layer thickness `[m]`"
    ΔZ::FT = ZS[1] - ZS[2]

    # Prognostic variables (not used for ∂y∂t)
    "Temperature `[K]`"
    t::FT = T₂₅(FT)

    # Prognostic variables (used for ∂y∂t)
    "Total stored energy per volume `[J m⁻³]`"
    e::FT = (CP * ρ + VC.Θ_SAT * CP_L(FT) * ρ_H₂O(FT)) * t
    "Soil water content"
    θ::FT = VC.Θ_SAT
    "Marginal increase in energy `[W m⁻²]`"
    ∂e∂t::FT = 0
    "Marginal increase in energy `[W m⁻²]`"
    ∂n∂t::Vector{FT} = zeros(FT,5)
    "Marginal increase in soil water content `[s⁻¹]`"
    ∂θ∂t::FT = 0

    # Diagnostic variables
    "Soil hydraulic conductance per area `[mol m⁻² s⁻¹ MPa⁻¹]`"
    k::FT = 0
    "Matric potential `[MPa]`"
    ψ::FT = 0

    # Cache variables
    "Combined specific heat capacity of soil `[J K⁻¹ kg⁻¹]`"
    _cp::FT = 0
    "Relative soil diffusive coefficient per area based on air fraction (distance accounted for already)"
    _kd::FT = 0
    "Combined soil thermal conductance `[W m⁻² K⁻¹]`"
    _λ_thermal::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jun-08: add Soil structure
#     2022-Jun-08: add constructor
#     2022-Jun-09: add fields: e_net_diffuse, e_net_direct
#     2022-Jun-10: add fields: r_net_lw, r_net_sw, ρ_lw
#     2022-Jun-14: add abstractized soil albedo
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jun-14: add field for soil color class
#     2022-Jul-13: move VC, Z, t, and θ to SoilLayer
#     2022-Jul-13: add field AREA, _k, _q, and _δψ
#     2022-Jul-13: add field _λ_thermal, _q_thermal, and _δt
#     2022-Jul-15: add field runoff
#     2023-Jun-16: remove fields DIM_*
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure for Soil

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Soil{FT<:AbstractFloat}
    # General information
    "Total area of the soil `[m²]`"
    AREA::FT = 500
    "Color class as in CLM"
    COLOR::Int = 1
    "Soil layers boundaries"
    ZS::Vector{FT} = FT[0,-0.1,-0.25,-0.5,-1,-3]

    # Embedded structures
    "Albedo related structure"
    ALBEDO::Union{BroadbandSoilAlbedo{FT}, HyperspectralSoilAlbedo{FT}}
    "Soil layers"
    LAYERS::Vector{SoilLayer{FT}} = SoilLayer{FT}[SoilLayer{FT}(VC = VanGenuchten{FT}("Loam"), ZS = ZS[_i:_i+1]) for _i in 1:length(ZS)-1]

    # Diagnostic variables
    "Surface runoff due to heavy precipitation during the time step `[mol m⁻²]`"
    runoff::FT = 0

    # Cache variables
    "Soil hydraulic conductance between layers per area `[mol m⁻² s⁻¹ MPa⁻¹]`"
    _k::Vector{FT}
    "Flux between layers per area `[mol m⁻² s⁻¹]`"
    _q::Vector{FT}
    "Thermal flux between layers per area `[mol m⁻² s⁻¹]`"
    _q_thermal::Vector{FT}
    "Soil temperature difference between layers `[MPa]`"
    _δt::Vector{FT}
    "Soil metric potential difference between layers `[MPa]`"
    _δψ::Vector{FT}
    "Soil thermal conductance between layers per area `[W m⁻² K⁻¹]`"
    _λ_thermal::Vector{FT}
end

Soil(config::SPACConfiguration{FT}; ground_area::Number = 500, soil_bounds::Vector{<:Number} = [0,-0.1,-0.25,-0.5,-1,-3]) where {FT} = (
    (; DIM_SOIL) = config;

    return Soil{FT}(
                AREA       = ground_area,
                ZS         = soil_bounds,
                ALBEDO     = HyperspectralSoilAlbedo(config),
                _k         = zeros(FT, DIM_SOIL - 1),
                _q         = zeros(FT, DIM_SOIL - 1),
                _q_thermal = zeros(FT, DIM_SOIL - 1),
                _δt        = zeros(FT, DIM_SOIL - 1),
                _δψ        = zeros(FT, DIM_SOIL - 1),
                _λ_thermal = zeros(FT, DIM_SOIL - 1),
    )
);
