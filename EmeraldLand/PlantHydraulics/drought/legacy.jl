#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-25: migrate inititialize_legacy! function to clear_legacy!
#     2022-May-26: add method for each hydraulic system like leaf, root, and stem
#     2022-May-26: add method for leaf (hydraulic system nested within)
#     2022-May-26: add method for SPAC system (hydraulic system nested within)
#     2022-Jun-29: rename SPAC to ML*SPAC to be more accurate
#     2022-Jun-30: add compatibility to Leaf
#     2022-Jul-08: deflate documentations
#     2024-Sep-03: use state.asap to check the xylem status (<= 0 means the xylem is dead)
#
#######################################################################################################################################################################################################
"""

    clear_legacy!(spac::BulkSPAC{FT}) where {FT}
    clear_legacy!(organ::Union{CanopyLayer{FT}, Leaf{FT}, Root{FT}, Stem{FT}}) where {FT}

Clear the legacy for hydraulic organ or system, given
- `spac` `BulkSPAC` type structure
- `organ` `Leaf`, `Root`, or `Stem` type structure
"""
function clear_legacy! end;

clear_legacy!(spac::BulkSPAC{FT}) where {FT} = (
    branches = spac.plant.branches;
    leaves = spac.plant.leaves;
    roots = spac.plant.roots;
    trunk = spac.plant.trunk;

    for root in roots
        clear_legacy!(root);
    end;

    clear_legacy!(trunk);

    for stem in branches
        clear_legacy!(stem);
    end;

    for leaf in leaves
        clear_legacy!(leaf);
    end;

    return nothing
);

clear_legacy!(organ::Union{CanopyLayer{FT}, Leaf{FT}, Root{FT}, Stem{FT}}) where {FT} = clear_legacy!(organ.xylem);

clear_legacy!(xylem::XylemHydraulics{FT}) where {FT} = (xylem.state.p_history .= 0; return nothing);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2024-Aug-31: migrate some features from plant_pressure_profile! to update_legacy!
#     2024-Sep-04: remove unnecessary methods (backward) for updating the legacy
#
#######################################################################################################################################################################################################
"""

    update_legacy!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update the legacy for hydraulic organ or system, given
- `config` `SPACConfiguration` type structure
- `spac` `BulkSPAC` type structure
"""
function update_legacy! end;

update_legacy!(
            config::SPACConfiguration{FT},
            spac::BulkSPAC{FT}) where {FT} = (
    for r in spac.plant.roots
        update_legacy!(config, r.xylem, r.energy.s_aux.t);
    end;
    update_legacy!(config, spac.plant.trunk.xylem, spac.plant.trunk.energy.s_aux.t);
    for s in spac.plant.branches
        update_legacy!(config, s.xylem, s.energy.s_aux.t);
    end;
    for l in spac.plant.leaves
        update_legacy!(config, l.xylem, l.energy.s_aux.t);
    end;

    return nothing
);

update_legacy!(
            config::SPACConfiguration{FT},
            xylem::XylemHydraulics{FT},
            t::FT) where {FT} = update_legacy!(config, xylem.state, xylem.auxil, t);

update_legacy!(
            config::SPACConfiguration{FT},
            x_state::XylemHydraulicsState{FT},
            x_aux::Union{XylemHydraulicsAuxilNSS{FT}, XylemHydraulicsAuxilSS{FT}},
            t::FT) where {FT} = (
    if x_state.asap <= 0 || !config.ENABLE_DROUGHT_LEGACY
        return nothing
    end;

    # update the pressure profile calculation only if xylem area > 0
    f_st = relative_surface_tension(t);
    N = length(x_state.p_history);
    for i in 1:N
        p_mem = x_state.p_history[i];
        p₂₅ = x_aux.pressure[i] / f_st;
        if p₂₅ < p_mem
            x_state.p_history[i] = p₂₅;
        end;
    end;

    return nothing
);
