#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-20: add function to create data from dict
#     2024-Feb-25: add t_aux! and s_aux! functions to update the axiliary variables
#
#######################################################################################################################################################################################################
"""

    spac_struct(gmdict::Dict, config::SPACConfiguration{FT}) where {FT}

Create a SPAC, given
- `gmdict` Dictionary of GriddingMachine data in a grid
- `config` Configurations for SPAC

"""
function spac_struct(gmdict::Dict, config::SPACConfiguration{FT}) where {FT}
    # read in canopy height
    z_canopy   = max(FT(0.1), gmdict["CANOPY_HEIGHT"]);
    Δz         = z_canopy / 20;
    air_bounds = Δz .* collect(0:21);

    # create a SPAC to work on
    spac = BulkSPAC(
                config;
                air_bounds = air_bounds,
                latitude = gmdict["LATITUDE"],
                longitude = gmdict["LONGITUDE"],
                soil_bounds = [0, -0.1, -0.35, -1, -3],
                plant_zs = [-2, z_canopy/2, z_canopy]);
    spac.soil_bulk.trait.color = gmdict["SOIL_COLOR"];

    # update soil type information per layer
    for i in eachindex(spac.soils)
        # TODO: add a line to parameterize K_MAX
        spac.soils[i].trait.vc.α = gmdict["SOIL_α"][i];
        spac.soils[i].trait.vc.N = gmdict["SOIL_N"][i];
        spac.soils[i].trait.vc.M = 1 - 1 / spac.soils[i].trait.vc.N;
        spac.soils[i].trait.vc.Θ_RES = gmdict["SOIL_ΘR"][i];
        spac.soils[i].trait.vc.Θ_SAT = gmdict["SOIL_ΘS"][i];
    end;

    # set hydraulic traits to very high so as to not triggering NaN (they do not impact result anyway)
    # for _organ in [spac.plant.leaves; spac.plant.branches; spac.plant.trunk; spac.plant.roots]
    #     _organ.xylem.trait.vc.B = 10;
    #     _organ.xylem.trait.vc.C = 1;
    # end;

    # update leaf mass per area and stomtal model
    for leaf in spac.plant.leaves
        leaf.bio.trait.lma = gmdict["LMA"];
    end;

    # TODO: add support to C4 photosynthesis
    # if gmdict["C3C4"] == "C4"
    #     error("C4 photosynthesis to be setted up");
    # end;

    # update the vcmax for C3 model
    prescribe_traits!(config, spac; vcmax = nanmean(gmdict["VCMAX25"]), vcmax_expo = 0.3);

    # initialize the spac
    initialize_states!(config, spac);
    initialize_spac!(config, spac);
    t_aux!(config, spac);
    s_aux!(config, spac);

    return spac
end;
