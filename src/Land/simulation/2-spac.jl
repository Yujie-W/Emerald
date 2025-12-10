"""

    site_spac(config::SPACConfig{FT}, gmd::Dict{String,Any}; lai_layer_strategy::Bool = true) where {FT}

Create a un-initialized SPAC using the data from a grid (CHL, VCMAX25, LAI, and CI are not prescribed as these changes with time), given
- `config` Configurations for SPAC
- `gmd` Dictionary of GriddingMachine data in a grid
- `lai_layer_strategy` Whether to use LAI-based air layer strategy (default: true)

"""
function site_spac(config::SPACConfig{FT}, gmd::Dict{String,Any}; lai_layer_strategy::Bool = true) where {FT}
    # compute air layer bounds based on maximum LAI
    zc = max(FT(0.05), gmd["CANOPY_HEIGHT"]);
    air_bounds = if lai_layer_strategy
        n_layer = max(10, Int(ceil(maximum(gmd["LAI"]) * 10)));
        collect(0:(2n_layer+2)) * zc / (2n_layer)
    else
        collect(0:21) * zc / 20
    end;

    #
    # TODO: add support to C4 photosynthesis
    #
    spac = BulkSPAC(
                config;
                air_bounds = air_bounds,
                elevation = gmd["ELEVATION"],
                latitude = gmd["LATITUDE"],
                longitude = gmd["LONGITUDE"],
                soil_bounds = [0, -0.1, -0.35, -1, -3],
                plant_zs = [-2, zc/2, zc]);
    # initialize soil color
    spac.soil_bulk.trait.color = gmd["SOIL_COLOR"];

    # set up LMA and infinite carbon pool (to avoid NSC pool)
    for i in eachindex(spac.plant.leaves)
        spac.plant.leaves[i].bio.trait.lma = gmd["LMA"];
    end;
    config.FEATURES.UNLIMITED_NSC_POOL ? spac.plant.pool.c_pool = Inf : nothing;

    # set up SAI
    prescribe_traits!(config, spac; lai = gmd["LAI"][1], sai = gmd["SAI"]);

    # update soil type information per layer
    for i in eachindex(spac.soils)
        # TODO: add a line to parameterize K_MAX
        # TODO: fix these later with better data source
        if !isnan(gmd["SOIL_α"][i]) && !isnan(gmd["SOIL_N"][i]) && !isnan(gmd["SOIL_ΘR"][i]) && !isnan(gmd["SOIL_ΘS"][i])
            spac.soils[i].trait.vc.α = gmd["SOIL_α"][i];
            spac.soils[i].trait.vc.N = gmd["SOIL_N"][i];
            spac.soils[i].trait.vc.M = 1 - 1 / spac.soils[i].trait.vc.N;
            spac.soils[i].trait.vc.Θ_RES = gmd["SOIL_ΘR"][i];
            spac.soils[i].trait.vc.Θ_SAT = gmd["SOIL_ΘS"][i];
            spac.soils[i].state.θ = spac.soils[i].trait.vc.Θ_SAT;
        end;
    end;

    # initialize the SPAC
    initialize_spac!(config, spac);

    return spac
end;
