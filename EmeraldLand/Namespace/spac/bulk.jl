# This file contains struct and constructor for the bulk SPAC system

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-17: add struct BulkSPAC
#     2023-Oct-17: add step to customize soil albedo algorithm from spac settings
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for monospecies tree SPAC system (with trunk and branches)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BulkSPAC{FT}
    "General information"
    info::SPACInfo{FT}

    "Soil bulk variables"
    soil_bulk::SoilBulk{FT}
    "Soil layers"
    soils::Vector{SoilLayer{FT}}

    "Air for each layer (more than canopy layer)"
    airs::Vector{AirLayer{FT}}
    "Meteorology information"
    meteo::Meteorology{FT}

    "Plant information"
    plant::Plant{FT}
    "Canopy used for radiation calculations"
    canopy::MultiLayerCanopy{FT}
end;

BulkSPAC(config::SPACConfiguration{FT};
         air_bounds::Vector{<:Number} = collect(0:0.5:13),
         basal_area::Number = 1,
         elevation::Number = 32,
         ground_area::Number = 500,
         latitude::Number = 33.173,
         longitude::Number = 115.4494,
         soil_bounds::Vector{<:Number} = [0,-0.1,-0.25,-0.5,-1,-3],
         plant_zs::Vector{<:Number} = [-1,6,12]
) where {FT} = (
    # set up the general information
    spac_info = SPACInfo{FT}(
                z_air  = air_bounds,
                z_soil = soil_bounds,
                elev   = elevation,
                lat    = latitude,
                lon    = longitude);

    # set up soil bulk parameters
    n_soil = length(soil_bounds) - 1;
    spac_sbulk = SoilBulk(config, n_soil);
    spac_sbulk.state.area = ground_area;

    # set up soil albedo algorithm
    if config.α_CLM && config.α_FITTING
        spac_sbulk.state.albedo = SoilAlbedoHyperspectralCLM();
    elseif config.α_CLM
        spac_sbulk.state.albedo = SoilAlbedoBroadbandCLM();
    elseif !config.α_CLM && config.α_FITTING
        spac_sbulk.state.albedo = SoilAlbedoHyperspectralCLIMA();
    else
        spac_sbulk.state.albedo = SoilAlbedoBroadbandCLIMA();
    end;

    # set up the soil layers (energy updated in initialize! function)
    soil_layers = SoilLayer{FT}[SoilLayer{FT}() for _ in 1:n_soil];
    for i in eachindex(soil_layers)
        soil_layers[i].state.zs = soil_bounds[i:i+1];
        soil_layers[i].auxil.z = (soil_bounds[i] + soil_bounds[i+1]) / 2;
        soil_layers[i].auxil.δz = (soil_bounds[i] - soil_bounds[i+1]);
    end;

    # set up the air layers
    n_air = length(air_bounds) - 1;
    air_layers = AirLayer{FT}[AirLayer{FT}() for _ in 1:n_air];
    for i in eachindex(air_layers)
        air_layers[i].state.zs = air_bounds[i:i+1];
        air_layers[i].auxil.z = (air_bounds[i] + air_bounds[i+1]) / 2;
        air_layers[i].auxil.δz = (air_bounds[i+1] - air_bounds[i]);
    end;

    # set up the meteorology
    spac_meteo = Meteorology{FT}(rad_sw = ShortwaveRadiation(config));

    # set up the plant (excluding canopy)
    # set up the roots
    mask_soil = findall(plant_zs[1] .< soil_bounds .< 0);
    n_root = length(mask_soil) + 1;
    ind_root = n_root > 1 ? [findfirst(plant_zs[1] .< soil_bounds .< 0) - 1; mask_soil] : mask_soil;
    roots = Root{FT}[Root(config) for _ in 1:n_root];
    for i in eachindex(roots)
        roots[i].xylem.state.area = basal_area / n_root;
        roots[i].xylem.state.Δh = 0 - max(plant_zs[1], soil_bounds[ind_root[i]+1]);
    end;

    # set up the trunk
    trunk = Stem(config);
    trunk.xylem.state.area = basal_area;
    trunk.xylem.state.Δh = plant_zs[2];

    # set up the branches
    mask_air = findall(plant_zs[2] .< air_bounds .< plant_zs[3]);
    n_layer = length(mask_air) + 1;
    ind_layer = n_layer > 1 ? [findfirst(plant_zs[2] .< air_bounds .< plant_zs[3])[1] - 1; mask_air] : mask_air;
    branches = Stem{FT}[Stem(config) for _ in 1:n_layer];
    for i in eachindex(branches)
        branches[i].xylem.state.area = basal_area / n_layer;
        branches[i].xylem.state.Δh = (min(plant_zs[3], air_bounds[ind_layer[i]+1]) - plant_zs[2]);
    end;

    # set up the leaves
    leaves = Leaf{FT}[Leaf(config) for i in 1:n_layer];

    # set up the plant
    plant = Plant{FT}(
                zs           = plant_zs,
                roots        = roots,
                roots_index  = ind_root,
                trunk        = trunk,
                branches     = branches,
                leaves       = leaves,
                leaves_index = ind_layer);

    # set up the canopy
    spac_canopy = MultiLayerCanopy(config, n_layer);

    return BulkSPAC{FT}(
                info      = spac_info,
                soil_bulk = spac_sbulk,
                soils     = soil_layers,
                airs      = air_layers,
                meteo     = spac_meteo,
                plant     = plant,
                canopy    = spac_canopy
    )
);
