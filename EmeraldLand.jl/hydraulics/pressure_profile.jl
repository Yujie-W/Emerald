
#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-27: migrate function to version v0.3
#     2022-Jul-12: compute e_crit for leaves
#     2022-Jul-12: compute β for leaves (only for empirical models)
#     2022-Jul-14: update root p_ups from SOIL
#     2022-Oct-20: use add SoilLayer to function variables, because of the removal of SH from RootHydraulics
#     2023-Sep-11: put option update to the SPAC configuration
#     2023-Sep-11: add config to the variable list
#     2023-Sep-11: rename function to critical_flow to xylem_pressure
# To do
#     TODO: add leaf extra-xylary vulnerability curve
#
#######################################################################################################################################################################################################
"""
This function is designed for the following purposes:
- Update organ pressure profile
- Update pressure profile for the entire SPAC

"""
function xylem_pressure_profile! end

# TODO: compute leaf e_crit
# TODO: update legacy
# TODO: make sure to not mixing with top soil that is meant for evaporation
# TODO: add soil ion concentration
"""

    xylem_pressure_profile!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Update xylem pressure profile (flow profile needs to be updated a priori), given
- `spac` `MultiLayerSPAC` type spac
- `drought_legacy` If true, update xylem cavitation legacy

"""
xylem_pressure_profile!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; KR_THRESHOLD) = config;
    (; BRANCHES, LEAVES, ROOTS, ROOTS_INDEX, SOIL, TRUNK) = spac;

    _nroots = length(ROOTS);

    # update water potential from SOIL
    for _i in eachindex(ROOTS_INDEX)
        _root = ROOTS[_i];
        _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];
        if _root._isconnected
            _root.HS.p_ups = soil_ψ_25(_slayer.VC, _slayer.θ) * relative_surface_tension(_slayer.t);
        else
            _root.HS.p_ups = xylem_pressure(_root.HS.VC, KR_THRESHOLD) * relative_surface_tension(_root.t);
        end;
    end;

    # update the profile in roots
    _p_mean::FT = 0;
    for _i in eachindex(ROOTS_INDEX)
        _root = ROOTS[_i];
        _slayer = SOIL.LAYERS[ROOTS_INDEX[_i]];
        xylem_pressure_profile!(config, _root, _slayer);
        _p_mean += (_root).HS.p_dos;
    end;
    _p_mean /= _nroots;

    # update the profile in trunk
    TRUNK.HS.p_ups = _p_mean;
    xylem_pressure_profile!(config, TRUNK);

    # update the profile in branch and leaf
    for _i in eachindex(BRANCHES)
        _stem = BRANCHES[_i];
        _leaf = LEAVES[_i];
        (_stem).HS.p_ups = TRUNK.HS.p_dos;
        xylem_pressure_profile!(config, _stem);
        (_leaf).HS.p_ups = (_stem).HS.p_dos;
        xylem_pressure_profile!(config, _leaf);
    end;

    # update the β factor for empirical models
    β_factor!(spac);

    return nothing
);
