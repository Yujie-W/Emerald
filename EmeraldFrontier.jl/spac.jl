
"""

    prepare_spac(dict::Dict; FT = Float64, c3c4::String = "C3", ePAR::Bool = true)

Create a SPAC, given
- `dict` Dictionary of GriddingMachine data in a grid
- `FT` Floating number type
- `c3c4` A string to indicate C3 or C4 photosynthesis
- `ePAR` Whether to use extended PAR in wl_set

"""
function prepare_spac(gm_dict::Dict)
    FT = gm_dict["FT"];
    # read in canopy height
    _z_canopy   = max(FT(0.1), gm_dict["CANOPY_HEIGHT"]);
    _Δz         = _z_canopy / 20;
    _air_bounds = _Δz .* collect(0:21);

    # create a SPAC to work on
    _spac = MultiLayerSPAC{FT}(
                DIM_AIR      = 25,
                DIM_LAYER    = 10,
                DIM_ROOT     = 4,
                LATITUDE     = gm_dict["LATITUDE"],
                LONGITUDE    = gm_dict["LONGITUDE"],
                LEAVES_INDEX = collect(11:20),
                ROOTS_INDEX  = collect(1:4),
                Z            = [-2, _z_canopy/2, _z_canopy],
                Z_AIR        = _air_bounds,
                SOIL         = Soil{FT}(DIM_SOIL = 4, COLOR = gm_dict["SOIL_COLOR"], ZS = [0, -0.1, -0.35, -1, -3]));

    # update soil type information per layer
    for _i in eachindex(_spac.SOIL.LAYERS)
        # TODO: add a line to parameterize K_MAX
        _spac.SOIL.LAYERS[_i].VC.α = gm_dict["SOIL_α"][_i];
        _spac.SOIL.LAYERS[_i].VC.N = gm_dict["SOIL_N"][_i];
        _spac.SOIL.LAYERS[_i].VC.M = 1 - 1 / _spac.SOIL.LAYERS[_i].VC.N;
        _spac.SOIL.LAYERS[_i].VC.Θ_RES = gm_dict["SOIL_Θr"][_i];
        _spac.SOIL.LAYERS[_i].VC.Θ_SAT = gm_dict["SOIL_Θs"][_i];
    end;

    # set hydraulic traits to very high so as to not triggering NaN (they do not impact result anyway)
    # for _organ in [_spac.LEAVES; _spac.BRANCHES; _spac.TRUNK; _spac.ROOTS]
    #     _organ.HS.VC.B = 10;
    #     _organ.HS.VC.C = 1;
    # end;

    # update leaf mass per area and stomtal model
    for _leaves in _spac.LEAVES
        _leaves.BIO.lma = gm_dict["LMA"];
    end;

    # add support to C4 photosynthesis
    if gm_dict["C3C4"] == "C4"
        @error "C4 photosynthesis to be setted up";
    end;

    # update the vcmax for C3 model
    update!(_spac; vcmax = nanmean(gm_dict["VCMAX25"]), vcmax_expo = 0.3);

    # initialize the spac
    initialize!(_spac);

    return _spac
end
