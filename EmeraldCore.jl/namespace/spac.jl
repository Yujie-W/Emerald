#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-11: add new struct for memory cache
#     2023-Mar-11: add new field temperature
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that store memory information

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SPACMemory{FT}
    "Chlorophyll content"
    chl::FT = -9999
    "Leaf area index"
    lai::FT = -9999
    "Temperature record for CLM T mean of 10 days"
    tem::Vector{FT} = ones(FT, 240) .* NaN
    "Vcmax25"
    vcm::FT = -9999
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-25: add abstract type for soil-plant-air continuum
#     2022-Jun-29: rename grass, palm, and tree SPAC to ML*SPAC
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractSPACSystem:
- [`MonoElementSPAC`](@ref)
- [`MultiLayerSPAC`](@ref)

"""
abstract type AbstractSPACSystem{FT} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: toy SPAC system
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-Jun-29: add AirLayer to SPAC
#     2022-Jul-14: add Meteorology to SPAC
#     2022-Mar-11: add MEMORY field
#     2023-Apr-13: add DIM_XYLEM to struct type
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for simplest SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MonoElementSPAC{FT,DIMS} <: AbstractSPACSystem{FT}
    # Embedded structures
    "Air conditions"
    AIR::AirLayer{FT} = AirLayer{FT}()
    "Leaf system"
    LEAF::Leaf{FT,DIMS} = Leaf{FT,DIMS}()
    "Memory cache"
    MEMORY::SPACMemory{FT} = SPACMemory{FT}()
    "Meteorology information"
    METEO::Meteorology{FT} = Meteorology{FT}()
    "Root system"
    ROOT::Root{FT,DIMS} = Root{FT,DIMS}()
    "Soil component"
    SOIL::Soil{FT,DIMS} = Soil{FT,DIMS}(ZS = FT[0, -1])
    "Stem system"
    STEM::Stem{FT,DIMS} = Stem{FT,DIMS}()

    # Cache variables
    "Relative hydraulic conductance"
    _krs::Vector{FT} = ones(FT, 4)
end

MonoElementSPAC(config::SPACConfiguration{FT,DIMS}) where {FT,DIMS} = MonoElementSPAC{FT,DIMS}();


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies tree
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-Jun-29: rename struct to MultiLayerSPAC, and use Leaves2D
#     2022-Jun-29: add CANOPY, Z, AIR, WLSET, LHA, ANGLES, SOIL, RAD_LW, RAD_SW, Φ_PHOTON to SPAC
#     2022-Jul-14: add Meteorology to SPAC
#     2022-Aug-30: remove LHA and WLSET
#     2023-Mar-11: add MEMORY and RAD_SW_REF fields
#     2023-Mar-28: add field _root_connection
#     2023-Apr-13: move Φ_PHOTON, RAD_SW_REF to SPACConfiguration
#     2023-Apr-13: move RAD_LW and RAD_SW to Meteorology
#     2023-Apr-13: add DIM_XYLEM to struct type
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies tree SPAC system (with trunk and branches)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MultiLayerSPAC{FT,DIMS} <: AbstractSPACSystem{FT}
    # Geometry information
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int} = collect(Int, 1:DIMS.DIM_CANOPY)
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int} = collect(Int, 1:DIMS.DIM_ROOT)
    "Depth and height information `[m]`"
    Z::Vector{FT} = FT[-1, 6, 12]
    "Air boundaries `[m]`"
    Z_AIR::Vector{FT} = collect(FT, 0:0.2:13)

    # Geographical information
    "Elevation"
    ELEVATION::FT = 32.0
    "Latitude"
    LATITUDE::FT = 33.173
    "Longitude"
    LONGITUDE::FT = 115.4494

    # Embedded structures
    "Air for each layer (more than canopy layer)"
    AIR::Vector{AirLayer{FT}} = AirLayer{FT}[AirLayer{FT}() for _i in 1:DIMS.DIM_AIR]
    "Sun sensor geometry"
    ANGLES::SunSensorGeometry{FT} = SunSensorGeometry{FT}()
    "Branch hydraulic system"
    BRANCHES::Vector{Stem{FT,DIMS}} = Stem{FT,DIMS}[Stem{FT,DIMS}() for _i in 1:DIMS.DIM_CANOPY]
    "Canopy used for radiation calculations"
    CANOPY::HyperspectralMLCanopy{FT,DIMS} = HyperspectralMLCanopy{FT,DIMS}()
    "Leaf per layer"
    LEAVES::Vector{Leaves2D{FT,DIMS}} = Leaves2D{FT,DIMS}[Leaves2D{FT,DIMS}() for _i in 1:DIMS.DIM_CANOPY]
    "Memory cache"
    MEMORY::SPACMemory{FT} = SPACMemory{FT}()
    "Meteorology information"
    METEO::Meteorology{FT,DIMS} = Meteorology{FT,DIMS}()
    "Root hydraulic system"
    ROOTS::Vector{Root{FT,DIMS}} = Root{FT,DIMS}[Root{FT,DIMS}() for _i in 1:DIMS.DIM_ROOT]
    "Soil component"
    SOIL::Soil{FT,DIMS} = Soil{FT,DIMS}()
    "Trunk hydraulic system"
    TRUNK::Stem{FT,DIMS} = Stem{FT,DIMS}()

    # Cache variables
    "Flow rate per root layer"
    _fs::Vector{FT} = zeros(FT, DIMS.DIM_ROOT)
    "Conductances for each root layer at given flow"
    _ks::Vector{FT} = zeros(FT, DIMS.DIM_ROOT)
    "Pressure for each root layer at given flow"
    _ps::Vector{FT} = zeros(FT, DIMS.DIM_ROOT)
    "Whether there is any root connected to soil"
    _root_connection::Bool = true
end

MultiLayerSPAC(config::SPACConfiguration{FT,DIMS}) where {FT,DIMS} = MultiLayerSPAC{FT,DIMS}();


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-13: add state struct to save
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for states of monospecies tree SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MultiLayerSPACState{FT,DIMS}
    # state variables
    "Shaded leaf stomatal conductance for all layers"
    gs_shaded::Vector{FT} = zeros(FT, DIMS.DIM_CANOPY)
    "Sunlit leaf stomatal conductance for all layers"
    gs_sunlit::Array{FT,3} = zeros(FT, DIMS.DIM_INCL, DIMS.DIM_AZI, DIMS.DIM_CANOPY)
    "Temperature record for CLM T mean of 10 days (based on CLM setting)"
    t_clm::Vector{FT} = ones(FT, DIMS.DIM_HISTORY)

    # variables to save
    "Gross primary productivity"
    gpp::FT = 0
    "TROPOMI SIF at 683 nm"
    tropomi_sif₆₈₃::FT = 0
    "TROPOMI SIF at 740 nm"
    tropomi_sif₇₄₀::FT = 0
end
