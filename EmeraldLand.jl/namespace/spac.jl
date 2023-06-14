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
Base.@kwdef mutable struct SPACMemory{FT<:AbstractFloat}
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
abstract type AbstractSPACSystem{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: toy SPAC system
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-Jun-29: add AirLayer to SPAC
#     2022-Jul-14: add Meteorology to SPAC
#     2022-Mar-11: add MEMORY field
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for simplest SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MonoElementSPAC{FT<:AbstractFloat} <: AbstractSPACSystem{FT}
    # Embedded structures
    "Air conditions"
    AIR::AirLayer{FT} = AirLayer{FT}()
    "Leaf system"
    LEAF::Leaf{FT} = Leaf{FT}()
    "Memory cache"
    MEMORY::SPACMemory{FT} = SPACMemory{FT}()
    "Meteorology information"
    METEO::Meteorology{FT} = Meteorology{FT}()
    "Root system"
    ROOT::Root{FT} = Root{FT}()
    "Soil component"
    SOIL::Soil{FT} = Soil{FT}(ZS = FT[0, -1])
    "Stem system"
    STEM::Stem{FT} = Stem{FT}()

    # Cache variables
    "Relative hydraulic conductance"
    _krs::Vector{FT} = ones(FT, 4)
end


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
#     2023-Apr-26: use a constructor rather than @kwdef
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies tree SPAC system (with trunk and branches)

# Fields

$(TYPEDFIELDS)

"""
mutable struct MultiLayerSPAC{FT<:AbstractFloat} <: AbstractSPACSystem{FT}
    # dimensions
    "Dimension of air layers"
    DIM_AIR::Int
    "Dimension of canopy layers"
    DIM_LAYER::Int
    "Dimension of root layers"
    DIM_ROOT::Int

    # Geometry information
    "Corresponding air layer per canopy layer"
    LEAVES_INDEX::Vector{Int}
    "Corresponding soil layer per root layer"
    ROOTS_INDEX::Vector{Int}
    "Depth and height information `[m]`"
    Z::Vector{FT}
    "Air boundaries `[m]`"
    Z_AIR::Vector{FT}

    # Geographical information
    "Elevation"
    ELEVATION::FT
    "Latitude"
    LATITUDE::FT
    "Longitude"
    LONGITUDE::FT

    # Embedded structures
    "Air for each layer (more than canopy layer)"
    AIR::Vector{AirLayer{FT}}
    "Sun sensor geometry"
    ANGLES::SunSensorGeometry{FT}
    "Branch hydraulic system"
    BRANCHES::Vector{Stem{FT}}
    "Canopy used for radiation calculations"
    CANOPY::HyperspectralMLCanopy{FT}
    "Leaf per layer"
    LEAVES::Vector{Leaves2D{FT}}
    "Memory cache"
    MEMORY::SPACMemory{FT}
    "Meteorology information"
    METEO::Meteorology{FT}
    "Root hydraulic system"
    ROOTS::Vector{Root{FT}}
    "Soil component"
    SOIL::Soil{FT}
    "Trunk hydraulic system"
    TRUNK::Stem{FT}

    # Cache variables
    "Flow rate per root layer"
    _fs::Vector{FT}
    "Conductances for each root layer at given flow"
    _ks::Vector{FT}
    "Pressure for each root layer at given flow"
    _ps::Vector{FT}
    "Whether there is any root connected to soil"
    _root_connection::Bool
end

MultiLayerSPAC{FT}(;
            air_bounds::Vector{<:Number} = collect(0:0.5:13),
            basal_area::Number = 1,
            elevation::Number = 32,
            ground_area::Number = 500,
            latitude::Number = 33.173,
            longitude::Number = 115.4494,
            soil_bounds::Vector{<:Number} = [0,-0.1,-0.25,-0.5,-1,-3],
            zs::Vector{<:Number} = [-1,6,12]
) where {FT<:AbstractFloat} = (
    _mask_air = findall(zs[2] .< air_bounds .< zs[3]);
    _mask_soil = findall(zs[1] .< soil_bounds .< 0);
    _dim_air = length(air_bounds) - 1;
    _dim_layer = length(_mask_air) + 1;
    _dim_root = length(_mask_soil) + 1;
    _dim_soil = length(soil_bounds) - 1;
    _ind_layer = _dim_layer > 1 ? [findfirst(zs[2] .< air_bounds .< zs[3])[1] - 1; _mask_air] : _mask_air;
    _ind_root = _dim_root > 1 ? [findfirst(zs[1] .< soil_bounds .< 0) - 1; _mask_soil] : _mask_soil;

    _air_layers = AirLayer{FT}[
        AirLayer{FT}(
            Z = (air_bounds[_i] + air_bounds[_i+1]) / 2,
            ΔZ = (air_bounds[_i+1] - air_bounds[_i])
        ) for _i in 1:_dim_air];
    _branches = Stem{FT}[
        Stem{FT}(
            HS = StemHydraulics{FT}(
                AREA = basal_area / _dim_layer,
                ΔH = (min(zs[3], air_bounds[_ind_layer[_i]+1]) - air_bounds[_ind_layer[_i]])
            )
        ) for _i in 1:_dim_layer];
    _roots = Root{FT}[
        Root{FT}(
            HS = RootHydraulics{FT}(
                AREA = basal_area / _dim_root,
                ΔH = (soil_bounds[_ind_root[_i]] - max(zs[1], soil_bounds[_ind_root[_i]+1]))
            )
        ) for _i in 1:_dim_root];

    return MultiLayerSPAC{FT}(
                _dim_air,                                                                   # DIM_AIR
                _dim_layer,                                                                 # DIM_LAYER
                _dim_root,                                                                  # DIM_ROOT
                _ind_layer,                                                                 # LEAVES_INDEX
                _ind_root,                                                                  # ROOTS_INDEX
                zs,                                                                         # Z
                air_bounds,                                                                 # Z_AIR
                elevation,                                                                  # ELEVATION
                latitude,                                                                   # LATITUDE
                longitude,                                                                  # LONGITUDE
                _air_layers,                                                                # AIR
                SunSensorGeometry{FT}(),                                                    # ANGLES
                _branches,                                                                  # BRANCHES
                HyperspectralMLCanopy{FT}(DIM_LAYER = _dim_layer),                          # CANOPY
                Leaves2D{FT}[Leaves2D{FT}() for _i in 1:_dim_layer],                        # LEAVES
                SPACMemory{FT}(),                                                           # MEMORY
                Meteorology{FT}(),                                                          # METEO
                _roots,                                                                     # ROOTS
                Soil{FT}(DIM_SOIL = _dim_soil, AREA = ground_area, ZS = soil_bounds),       # SOIL
                Stem{FT}(HS = StemHydraulics{FT}(AREA = basal_area, ΔH = zs[2] - zs[1])),   # TRUNK
                zeros(FT, _dim_root),                                                       # _fs
                zeros(FT, _dim_root),                                                       # _ks
                zeros(FT, _dim_root),                                                       # _ps
                true,                                                                       # _root_connection
    )
);


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
Base.@kwdef mutable struct MultiLayerSPACState{FT}
    # state variables
    "Shaded leaf stomatal conductance for all layers"
    gs_shaded::Vector{FT}
    "Sunlit leaf stomatal conductance for all layers"
    gs_sunlit::Array{FT,3}
    "Temperature record for CLM T mean of 10 days (based on CLM setting)"
    t_clm::Vector{FT}

    # variables to save
    "Soil moisture tuning factor β"
    beta::FT = 0
    "Gross primary productivity"
    gpp::FT = 0
    "MODIS EVI"
    modis_evi::FT = 0
    "MODIS NDVI"
    modis_ndvi::FT = 0
    "MODIS NIRv"
    modis_nirv::FT = 0
    "OCO SIF at 759 nm"
    oco_sif₇₅₉::FT = 0
    "OCO SIF at 770 nm"
    oco_sif₇₇₀::FT = 0
    "PAR"
    par::FT = 0
    "PPAR"
    ppar::FT = 0
    "Transpiration"
    transpiration::FT = 0
    "TROPOMI SIF at 683 nm"
    tropomi_sif₆₈₃::FT = 0
    "TROPOMI SIF at 740 nm"
    tropomi_sif₇₄₀::FT = 0
end

MultiLayerSPACState{FT}(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; DIM_LAYER, LEAVES) = spac;

    _gs_sunlit = zeros(FT, LEAVES[1].DIM_INCL, LEAVES[1].DIM_AZI, DIM_LAYER);
    for _i in 1:DIM_LAYER
        _gs_sunlit[:,:,_i] .= LEAVES[_i].g_H₂O_s_sunlit;
    end;

    return MultiLayerSPACState{FT}(
                gs_shaded = [_leaves.g_H₂O_s_shaded for _leaves in LEAVES],
                gs_sunlit = _gs_sunlit,
                t_clm = deepcopy(spac.MEMORY.tem),
    )
);
