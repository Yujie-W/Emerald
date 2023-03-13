#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the CACHE_SPAC
#
#######################################################################################################################################################################################################
"""

    initialize_cache!(FT)

Initialize the global parameter `CACHE_SPAC`, given
- `FT` Floating type

"""
function initialize_cache!(FT)
    global CACHE_SPAC;

    # create a SPAC to work on
    _z_canopy = FT(10);
    CACHE_SPAC = MonoMLTreeSPAC{FT}(
                DIM_AIR      = 25,
                DIM_LAYER    = 10,
                DIM_ROOT     = 4,
                LATITUDE     = 0,
                LONGITUDE    = 0,
                LEAVES_INDEX = collect(11:20),
                ROOTS_INDEX  = collect(1:4),
                Z            = [-2, _z_canopy/2, _z_canopy],
                Z_AIR        = collect(0:21) * _z_canopy / 20,
                SOIL         = Soil{FT}(DIM_SOIL = 4, ZS = [0, -0.1, -0.35, -1, -3]));

    # set hydraulic traits to very high so as to not triggering NaN (they do not impact result anyway)
    for _organ in [CACHE_SPAC.LEAVES; CACHE_SPAC.BRANCHES; CACHE_SPAC.TRUNK; CACHE_SPAC.ROOTS]
        _organ.HS.VC.B = 3;
        _organ.HS.VC.C = 1;
    end;

    # update leaf mass per area and stomtal model
    @inline linear_p_soil(x) = min(1, max(eps(FT), 1 + x / 5));
    _bt = BetaFunction{FT}(FUNC = linear_p_soil, PARAM_X = BetaParameterPsoil(), PARAM_Y = BetaParameterG1());
    for _leaves in CACHE_SPAC.LEAVES
        _leaves.SM = MedlynSM{FT}(G0 = 0.005, β = _bt);
    end;

    # initialize the spac
    initialize!(CACHE_SPAC);

    return nothing
end;
