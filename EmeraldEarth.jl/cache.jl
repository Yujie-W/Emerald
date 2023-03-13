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


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the CACHE_SPAC
#
#######################################################################################################################################################################################################
"""

    synchronize_cache!(gm_params::Nothing)
    synchronize_cache!(gm_params::Dict{String,Any})

Synchronize SPAC parameters from,
- `gm_params` Dict for GriddingMachine parameters, or nothing (if not land)

"""
function synchronize_cache! end

synchronize_cache!(gm_params::Nothing) = (return nothing);

synchronize_cache!(gm_params::Dict{String,Any}) = (
    FT = gm_params["FT"];
    _z_canopy = max(FT(0.1), gm_params["CANOPY_HEIGHT"]);

    #
    # TODO: update canopy height for plant hydraulic system based on _z_canopy
    #

    # update the values in the CACHE_SPAC, use .= for arrays
    global CACHE_SPAC;
    CACHE_SPAC.LATITUDE = gm_params["LATITUDE"];
    CACHE_SPAC.LONGITUDE = gm_params["LONGITUDE"];
    CACHE_SPAC.Z .= [-2, _z_canopy/2, _z_canopy];
    CACHE_SPAC.Z_AIR .= collect(0:21) * _z_canopy / 20;
    CACHE_SPAC.SOIL.COLOR = gm_params["SOIL_COLOR"];

    # update soil type information per layer
    for _i in eachindex(CACHE_SPAC.SOIL.LAYERS)
        # TODO: add a line to parameterize K_MAX
        CACHE_SPAC.SOIL.LAYERS[_i].VC.α = gm_params["SOIL_α"][_i];
        CACHE_SPAC.SOIL.LAYERS[_i].VC.N = gm_params["SOIL_N"][_i];
        CACHE_SPAC.SOIL.LAYERS[_i].VC.M = 1 - 1 / CACHE_SPAC.SOIL.LAYERS[_i].VC.N;
        CACHE_SPAC.SOIL.LAYERS[_i].VC.Θ_RES = gm_params["SOIL_ΘR"][_i];
        CACHE_SPAC.SOIL.LAYERS[_i].VC.Θ_SAT = gm_params["SOIL_ΘS"][_i];
    end;

    # update leaf mass per area and stomtal model
    for _leaves in CACHE_SPAC.LEAVES
        _leaves.BIO.lma = gm_params["LMA"];
        _leaves.SM.G1 = gm_params["MEDLYN_G1"];
    end;

    # sync the environmental conditions per layer for CO₂ concentration
    if !isnothing(gm_params["CO2"])
        for _alayer in CACHE_SPAC.AIR
            update!(_alayer; f_CO₂ = gm_params["CO2"]);
        end;
    end;

    return nothing
);
