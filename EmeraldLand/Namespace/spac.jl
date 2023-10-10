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
# Changes to this struct
# General
#     2022-May-25: SPAC system for monospecies tree
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-Jun-29: rename struct to MultiLayerSPAC, and use Leaf
#     2022-Jun-29: add CANOPY, Z, AIRS, WLSET, LHA, ANGLES, SOIL, RAD_LW, RAD_SW, Φ_PHOTON to SPAC
#     2022-Jul-14: add Meteorology to SPAC
#     2022-Aug-30: remove LHA and WLSET
#     2023-Mar-11: add MEMORY and RAD_SW_REF fields
#     2023-Mar-28: add field _root_connection
#     2023-Apr-13: move Φ_PHOTON, RAD_SW_REF to SPACConfiguration
#     2023-Apr-13: move RAD_LW and RAD_SW to Meteorology
#     2023-Apr-26: use a constructor rather than @kwdef
#     2023-Jun-16: remove fields DIM_*
#     2023-Aug-27: make corrections over the delta height of the xylem
#     2023-Sep-28: add field JUNCTION
#     2023-Oct-05: add fields SOIL_BULK and SOILS
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies tree SPAC system (with trunk and branches)

# Fields

$(TYPEDFIELDS)

"""
mutable struct MultiLayerSPAC{FT}
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
    AIRS::Vector{AirLayer{FT}}
    "Sun sensor geometry"
    ANGLES::SunSensorGeometry{FT}
    "Branch hydraulic system"
    BRANCHES::Vector{Stem{FT}}
    "Canopy used for radiation calculations"
    CANOPY::HyperspectralMLCanopy{FT}
    "Leaf per layer"
    LEAVES::Vector{Leaf{FT}}
    "Memory cache"
    MEMORY::SPACMemory{FT}
    "Meteorology information"
    METEO::Meteorology{FT}
    "Root hydraulic system"
    ROOTS::Vector{Root{FT}}
    "Soil bulk variables"
    SOIL_BULK::SoilBulk{FT}
    "Soil layers"
    SOILS::Vector{SoilLayer{FT}}
    "Trunk hydraulic system"
    TRUNK::Stem{FT}
    "Root-trunk junction capacitor used for roots flow calculations"
    JUNCTION::JunctionCapacitor{FT}

    # Cache variables
    "Whether there is any root connected to soil"
    _root_connection::Bool
end

MultiLayerSPAC(
            config::SPACConfiguration{FT};
            air_bounds::Vector{<:Number} = collect(0:0.5:13),
            basal_area::Number = 1,
            elevation::Number = 32,
            ground_area::Number = 500,
            latitude::Number = 33.173,
            longitude::Number = 115.4494,
            soil_bounds::Vector{<:Number} = [0,-0.1,-0.25,-0.5,-1,-3],
            zs::Vector{<:Number} = [-1,6,12]
) where {FT} = (
    mask_air = findall(zs[2] .< air_bounds .< zs[3]);
    mask_soil = findall(zs[1] .< soil_bounds .< 0);
    config.DIM_AIR = length(air_bounds) - 1;
    config.DIM_LAYER = length(mask_air) + 1;
    config.DIM_ROOT = length(mask_soil) + 1;
    config.DIM_SOIL = length(soil_bounds) - 1;
    ind_layer = config.DIM_LAYER > 1 ? [findfirst(zs[2] .< air_bounds .< zs[3])[1] - 1; mask_air] : mask_air;
    ind_root = config.DIM_ROOT > 1 ? [findfirst(zs[1] .< soil_bounds .< 0) - 1; mask_soil] : mask_soil;

    # set up the soil layers (energy updated in initialize! function)
    soil_layers = SoilLayer{FT}[SoilLayer{FT}() for _ in 1:config.DIM_SOIL];
    for i in eachindex(soil_layers)
        soil_layers[i].state.zs = soil_bounds[i:i+1];
        soil_layers[i].auxil.z = (soil_bounds[i] + soil_bounds[i+1]) / 2;
        soil_layers[i].auxil.δz = (soil_bounds[i] - soil_bounds[i+1]);
    end;

    # set up the roots
    roots = Root{FT}[Root(config) for _ in 1:config.DIM_ROOT];
    for i in eachindex(roots)
        roots[i].xylem.state.area = basal_area / config.DIM_ROOT;
        roots[i].xylem.state.Δh = 0 - max(zs[1], soil_bounds[ind_root[i]+1]);
    end;

    # set up the trunk
    trunk = Stem(config);
    trunk.xylem.state.area = basal_area;
    trunk.xylem.state.Δh = zs[2];

    # set up the branches
    branches = Stem{FT}[Stem(config) for _ in 1:config.DIM_LAYER];
    for i in eachindex(branches)
        branches[i].xylem.state.area = basal_area / config.DIM_LAYER;
        branches[i].xylem.state.Δh = (min(zs[3], air_bounds[ind_layer[i]+1]) - zs[2]);
    end;

    # set up the air layers
    air_layers = AirLayer{FT}[AirLayer{FT}() for _ in 1:config.DIM_AIR];
    for i in eachindex(air_layers)
        air_layers[i].state.zs = air_bounds[i:i+1];
        air_layers[i].auxil.z = (air_bounds[i] + air_bounds[i+1]) / 2;
        air_layers[i].auxil.δz = (air_bounds[i+1] - air_bounds[i]);
    end;

    return MultiLayerSPAC{FT}(
                ind_layer,                                              # LEAVES_INDEX
                ind_root,                                               # ROOTS_INDEX
                zs,                                                     # Z
                air_bounds,                                             # Z_AIR
                elevation,                                              # ELEVATION
                latitude,                                               # LATITUDE
                longitude,                                              # LONGITUDE
                air_layers,                                             # AIRS
                SunSensorGeometry{FT}(),                                # ANGLES
                branches,                                               # BRANCHES
                HyperspectralMLCanopy(config),                          # CANOPY
                Leaf{FT}[Leaf(config) for i in 1:config.DIM_LAYER],     # LEAVES
                SPACMemory{FT}(),                                       # MEMORY
                Meteorology{FT}(rad_sw = ShortwaveRadiation(config)),   # METEO
                roots,                                                  # ROOTS
                SoilBulk(config),                                       # SOIL_BULK
                soil_layers,                                            # SOILS
                trunk,                                                  # TRUNK
                JunctionCapacitor{FT}(),                                # JUNCTION
                true,                                                   # _root_connection
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
    "Total chloroplast fluorescence"
    csif::FT = 0
    "Total ETR"
    etr::FT = 0
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

MultiLayerSPACState{FT}(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; LEAVES) = spac;

    _gs_sunlit = zeros(FT, size(LEAVES[1].g_H₂O_s_sunlit,1), size(LEAVES[1].g_H₂O_s_sunlit,2), length(LEAVES));
    for i in eachindex(LEAVES)
        _gs_sunlit[:,:,i] .= LEAVES[i].g_H₂O_s_sunlit;
    end;

    return MultiLayerSPACState{FT}(
                gs_shaded = [_leaves.g_H₂O_s_shaded for _leaves in LEAVES],
                gs_sunlit = _gs_sunlit,
                t_clm = deepcopy(spac.MEMORY.tem),
    )
);
