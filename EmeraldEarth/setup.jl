#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the CACHE_SPAC
#     2023-Mar-13: initialize CACHE_STATE at the same time
#     2023-Mar-13: initialize CACHE_CONFIG at the same time
#     2023-Jun-15: make sure prescribed swc does not exceed the limits
#     2024-Feb-22: remove state and auxil from spac struct
#     2024-Feb-23: rename to setup_cache!
# Bug fixes
#     2023-Aug-26: make sure sza < 89 when total radiation is higher than 10 W m⁻²
#
#######################################################################################################################################################################################################
"""

    setup_cache!(FT::DataType = Float64)

Initialize the global parameter `CACHE_SPAC` (in all threads after loading workers), given
- `FT` Floating type (default is Float64)

"""
function setup_cache!(FT::DataType = Float64)
    global CACHE_CONFIG, CACHE_SPAC, CACHE_STATE;

    # create a SPAC to work on
    z_canopy = FT(10);
    CACHE_CONFIG = SPACConfiguration{FT}();
    CACHE_SPAC = BulkSPAC(
                CACHE_CONFIG;
                air_bounds = collect(0:21) * z_canopy / 20,
                latitude = 0,
                longitude = 0,
                soil_bounds = [0, -0.1, -0.35, -1, -3],
                plant_zs = [-2, z_canopy/2, z_canopy]);

    # update leaf mass per area and stomtal model
    @inline linear_p_soil(x) = min(1, max(eps(FT), 1 + x / 5));
    bt = BetaFunction{FT}(FUNC = linear_p_soil, PARAM_X = BetaParameterPsoil(), PARAM_Y = BetaParameterG1());
    for leaf in CACHE_SPAC.plant.leaves
        leaf.flux.state.stomatal_model = MedlynSM{FT}(G0 = 0.005, β = bt);
    end;

    # initialize the spac with non-saturated soil
    prescribe_soil!(CACHE_SPAC; swcs = Tuple(max(soil.state.vc.Θ_SAT - 0.02, (soil.state.vc.Θ_SAT + soil.state.vc.Θ_RES) / 2) for soil in CACHE_SPAC.soils));
    initialize!(CACHE_CONFIG, CACHE_SPAC);

    # create a state struct based on the spac
    CACHE_STATE = BulkSPACStates(CACHE_SPAC);

    return nothing
end;
