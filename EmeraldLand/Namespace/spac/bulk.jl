# This file contains struct and constructor for the bulk SPAC system

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-17: add struct BulkSPAC
#     2023-Oct-17: add step to customize soil albedo algorithm from spac settings
#     2024-Feb-27: customize leaf area from leaf area index per layer in the constructor function
#     2024-Jul-24: add field cache to store cache information for speeding up the model
#     2024-Aug-05: add Y = 1 - beta ^ (100 * -z) to calculate the root area distribution
#     2024-Aug-29: set crown radius and height change as xylem length
#     2024-Aug-30: initialize carbon pool for the whole plant when using the constructor
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

    # Cache structure
    "Cache structure for the SPAC model"
    cache::SPACCache{FT}
end;

BulkSPAC(config::SPACConfiguration{FT};
         air_bounds::Vector{<:Number} = collect(0:0.5:13),
         basal_area::Number = 0.15,
         elevation::Number = 32,
         ground_area::Number = 80,
         latitude::Number = 33.173,
         longitude::Number = 115.4494,
         max_lai::Number = 3,
         root_beta::Number = 0.961,
         soil_bounds::Vector{<:Number} = [0,-0.1,-0.25,-0.5,-1,-3],
         plant_zs::Vector{<:Number} = [-1,6,12],
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
    spac_sbulk.trait.area = ground_area;

    # set up the soil layers (energy updated in initialize! function)
    soil_layers = SoilLayer{FT}[SoilLayer{FT}() for _ in 1:n_soil];
    for i in eachindex(soil_layers)
        soil_layers[i].trait.zs = soil_bounds[i:i+1];
    end;

    # set up the air layers
    n_air = length(air_bounds) - 1;
    air_layers = AirLayer{FT}[AirLayer{FT}() for _ in 1:n_air];
    for i in eachindex(air_layers)
        air_layers[i].trait.zs = air_bounds[i:i+1];
    end;

    # set up the meteorology
    spac_meteo = Meteorology{FT}(rad_sw = ShortwaveRadiation(config));

    # set up the plant (excluding canopy)
    crown_radius = sqrt(ground_area / π);

    # set up the roots
    mask_soil = findall(plant_zs[1] .< soil_bounds .< 0);
    n_root = length(mask_soil) + 1;
    ind_root = n_root > 1 ? [findfirst(plant_zs[1] .< soil_bounds .< 0) - 1; mask_soil] : mask_soil;
    roots = Root{FT}[Root(config) for _ in 1:n_root];
    normi = 1 - root_beta ^ (100 * -soil_bounds[n_root+1]);
    for i in eachindex(roots)
        roots[i].xylem.trait.area = basal_area * (root_beta ^ (100 * -soil_bounds[ind_root[i]]) - root_beta ^ (100 * -soil_bounds[ind_root[i]+1])) / normi;
        roots[i].xylem.state.asap = roots[i].xylem.trait.area;
        roots[i].xylem.trait.Δh = 0 - max(plant_zs[1], soil_bounds[ind_root[i]+1]);
        roots[i].xylem.trait.l = crown_radius + roots[i].xylem.trait.Δh;
    end;

    # set up the trunk
    trunk = Stem(config);
    trunk.xylem.trait.area = basal_area;
    trunk.xylem.state.asap = trunk.xylem.trait.area;
    trunk.xylem.trait.Δh = plant_zs[2];
    trunk.xylem.trait.l = trunk.xylem.trait.Δh;

    # set up the branches
    mask_air = findall(plant_zs[2] .< air_bounds .< plant_zs[3]);
    n_layer = length(mask_air) + 1;
    ind_layer = n_layer > 1 ? [findfirst(plant_zs[2] .< air_bounds .< plant_zs[3])[1] - 1; mask_air] : mask_air;
    branches = Stem{FT}[Stem(config) for _ in 1:n_layer];
    for i in eachindex(branches)
        branches[i].xylem.trait.area = basal_area / n_layer;
        branches[i].xylem.state.asap = branches[i].xylem.trait.area;
        branches[i].xylem.trait.Δh = (min(plant_zs[3], air_bounds[ind_layer[i]+1]) - plant_zs[2]);
        branches[i].xylem.trait.l = crown_radius + branches[i].xylem.trait.Δh;
    end;

    # set up the canopy
    spac_canopy = MultiLayerCanopy(config, n_layer);

    # set up the canopy layers
    leaves = CanopyLayer{FT}[CanopyLayer(config) for i in 1:n_layer];
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaves[ilf].xylem.trait.area = spac_sbulk.trait.area * spac_canopy.structure.trait.δlai[irt];
        leaves[ilf].xylem.state.asap = leaves[ilf].xylem.trait.area;
    end;

    # set up the carbon pools
    lai_pool = ground_area * max_lai * leaves[1].bio.trait.lma * 10000 / 30;
    lai_pool_max = 2.5 * lai_pool;
    c_pool = CarbonPoolWholePlant{FT}(lai_pool, lai_pool_max);

    # set up the plant
    plant = Plant{FT}(
                zs           = plant_zs,
                z_beta       = root_beta,
                roots        = roots,
                roots_index  = ind_root,
                trunk        = trunk,
                branches     = branches,
                leaves       = leaves,
                leaves_index = ind_layer,
                pool         = c_pool,
                memory       = PlantMemory(config));

    # set up the cache
    cache = SPACCache{FT}(
                config.DIM_AZI,
                config.DIM_INCL,
                n_layer,
                config.DIM_PPAR_BINS,
                length(config.SPECTRA.Λ_SIF),
                length(config.SPECTRA.Λ_SIFE),
                length(config.SPECTRA.Λ));

    return BulkSPAC{FT}(
                info      = spac_info,
                soil_bulk = spac_sbulk,
                soils     = soil_layers,
                airs      = air_layers,
                meteo     = spac_meteo,
                plant     = plant,
                canopy    = spac_canopy,
                cache     = cache
    )
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-2: add struct BulkSPACStates
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to save bulk SPAC states (collections of states)

# Fields

$(TYPEDFIELDS)

"""
mutable struct BulkSPACStates{FT}
    "Soil layers"
    soils::Vector{SoilLayerState{FT}}

    "Air for each layer (more than canopy layer)"
    airs::Vector{AirLayerState{FT}}

    "Plant information"
    plant::PlantStates{FT}
    "Canopy used for radiation calculations"
    canopy::MultiLayerCanopyStates{FT}
end;

BulkSPACStates(spac::BulkSPAC{FT}) where {FT} = BulkSPACStates{FT}(
            [deepcopy(soil.state) for soil in spac.soils],
            [deepcopy(air.state) for air in spac.airs],
            PlantStates(spac.plant),
            MultiLayerCanopyStates(spac.canopy)
);

sync_state!(spac::BulkSPAC{FT}, states::BulkSPACStates{FT}) where {FT} = (
    for i in eachindex(spac.soils)
        sync_state!(spac.soils[i].state, states.soils[i]);
    end;
    for i in eachindex(spac.airs)
        sync_state!(spac.airs[i].state, states.airs[i]);
    end;
    sync_state!(spac.plant, states.plant);
    sync_state!(spac.canopy, states.canopy);

    return nothing
);

sync_state!(states::BulkSPACStates{FT}, spac::BulkSPAC{FT}) where {FT} = (
    for i in eachindex(states.soils)
        sync_state!(states.soils[i], spac.soils[i].state);
    end;
    for i in eachindex(states.airs)
        sync_state!(states.airs[i], spac.airs[i].state);
    end;
    sync_state!(states.plant, spac.plant);
    sync_state!(states.canopy, spac.canopy);

    return nothing
);
