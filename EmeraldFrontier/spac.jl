#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-20: add function to create data from dict
#
#######################################################################################################################################################################################################
"""

    spac(gmdict::Dict, config::SPACConfiguration{FT}) where {FT}

Create a SPAC, given
- `gmdict` Dictionary of GriddingMachine data in a grid
- `config` Configurations for SPAC

"""
function spac(gmdict::Dict, config::SPACConfiguration{FT}) where {FT}
    # read in canopy height
    _z_canopy   = max(FT(0.1), gmdict["CANOPY_HEIGHT"]);
    _Δz         = _z_canopy / 20;
    _air_bounds = _Δz .* collect(0:21);

    # create a SPAC to work on
    _spac = MultiLayerSPAC(
                config;
                air_bounds = _air_bounds,
                latitude = gmdict["LATITUDE"],
                longitude = gmdict["LONGITUDE"],
                soil_bounds = [0, -0.1, -0.35, -1, -3],
                zs = [-2, _z_canopy/2, _z_canopy]);
    _spac.SOIL_BULK.state.color = gmdict["SOIL_COLOR"];

    # update soil type information per layer
    for i in eachindex(_spac.SOILS)
        # TODO: add a line to parameterize K_MAX
        _spac.SOILS[i].state.vc.α = gmdict["SOIL_α"][i];
        _spac.SOILS[i].state.vc.N = gmdict["SOIL_N"][i];
        _spac.SOILS[i].state.vc.M = 1 - 1 / _spac.SOILS[i].state.vc.N;
        _spac.SOILS[i].state.vc.Θ_RES = gmdict["SOIL_ΘR"][i];
        _spac.SOILS[i].state.vc.Θ_SAT = gmdict["SOIL_ΘS"][i];
    end;

    # set hydraulic traits to very high so as to not triggering NaN (they do not impact result anyway)
    # for _organ in [_spac.LEAVES; _spac.BRANCHES; _spac.TRUNK; _spac.ROOTS]
    #     _organ.xylem.state.vc.B = 10;
    #     _organ.xylem.state.vc.C = 1;
    # end;

    # update leaf mass per area and stomtal model
    for _leaves in _spac.LEAVES
        _leaves.bio.state.lma = gmdict["LMA"];
    end;

    # add support to C4 photosynthesis
    if gmdict["C3C4"] == "C4"
        error("C4 photosynthesis to be setted up");
    end;

    # update the vcmax for C3 model
    update!(config, _spac; vcmax = nanmean(gmdict["VCMAX25"]), vcmax_expo = 0.3);

    # initialize the spac
    initialize!(config, _spac);

    return _spac
end
