# This file contains function to test what happens if the xylem hydraulic system is in a recovery manner (not realistic, but has to be done this before I use segmented xylem VC)

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2024-Aug-30: add function to determine whether to grow new xylem or recover the old one (function-wise, not real recovery)
#
#######################################################################################################################################################################################################
"""

    recovery_or_growth(organ::Union{Root{FT}, Stem{FT}}, c_mol::FT) where {FT}

Determine whether to grow new xylem or recover the old one (function-wise, not real recovery), given
- `organ` `Root` or `Stem` type structure
- `c_mol` Carbon investment `[mol]`

"""
function recovery_or_growth end;

recovery_or_growth(organ::Union{Root{FT}, Stem{FT}}, c_mol::FT) where {FT} = recovery_or_growth(organ.xylem, c_mol, organ.energy.s_aux.t);

recovery_or_growth(xylem::XylemHydraulics{FT}, c_mol::FT, t::FT) where {FT} = (
    f_st = relative_surface_tension(t);
    N = length(xylem.state.p_history);

    # compute the potential area from the given carbon investment
    delta_a = c_mol / (xylem.trait.l * xylem.trait.ρ * 1000 / 30);

    # compute the conductance before recovery
    Σr = 0;
    for i in 1:N
        k_i = xylem.state.asap * relative_xylem_k(xylem.trait.vc, xylem.state.p_history[i]) * xylem.trait.k_max * N;
        Σr += 1 / k_i;
    end;
    k_0 = 1 / Σr;

    # compute the conductance after recovery
    Σr = 0;
    for i in 1:N
        relative_k_his = relative_xylem_k(xylem.trait.vc, xylem.state.p_history[i]) + delta_a / xylem.state.asap;
        relative_k_cur = relative_xylem_k(xylem.trait.vc, xylem.auxil.pressure[i] / f_st);
        k_i = xylem.state.asap * min(relative_k_cur, relative_k_his) * xylem.trait.k_max * N;
        Σr += 1 / k_i;
    end;
    k_1 = 1 / Σr;

    # compute the conductance after new growth
    k_2 = k_0 * (1 + delta_a / xylem.state.asap);

    return k_1 > k_2
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2024-Aug-30: add function to recover the xylem hydraulic system (area added to total area, pressure history updated)
#     2024-Sep-04: add method to recover leaf xylem hydraulic system from new growth
#     2024-Oct-20: reconnect leaf when recover the xylem hydraulic system for leaf
#
#######################################################################################################################################################################################################
"""

    xylem_recovery!(organ::Union{Root{FT}, Stem{FT}}, c_mol::FT) where {FT}

Recover the xylem hydraulic system (area added to total area, pressure history updated), given
- `organ` `Root` or `Stem` type structure
- `c_mol` Carbon investment `[mol]`

"""
function xylem_recovery! end;

xylem_recovery!(organ::Union{Root{FT}, Stem{FT}}, c_mol::FT) where {FT} = xylem_recovery!(organ.xylem, c_mol);

# for stem and root only
xylem_recovery!(xylem::XylemHydraulics{FT}, c_mol::FT) where {FT} = (
    N = length(xylem.state.p_history);

    # compute the potential area from the given carbon investment
    delta_a = c_mol / (xylem.trait.l * xylem.trait.ρ * 1000 / 30);

    # compute the conductance after recovery
    for i in 1:N
        k_rel = relative_xylem_k(xylem.trait.vc, xylem.state.p_history[i]) + delta_a / xylem.state.asap;
        p_his = xylem_pressure(xylem.trait.vc, min(1, k_rel));
        xylem.state.p_history[i] = p_his;
    end;

    # add the new xylem area to the total area
    xylem.trait.area += delta_a;

    return nothing
);

# for leaf only
xylem_recovery!(xylem::XylemHydraulics{FT}, lai_0::FT, lai_diff::FT) where {FT} = (
    @assert lai_diff > 0 "The difference in LAI must be positive";

    # compute the conductance after recovery
    N = length(xylem.state.p_history);
    for i in 1:N
        k_rel_0 = relative_xylem_k(xylem.trait.vc, xylem.state.p_history[i]);
        k_rel_1 = (k_rel_0 * lai_0 + lai_diff) / (lai_0 + lai_diff);
        p_his = xylem_pressure(xylem.trait.vc, min(1, k_rel_1));
        xylem.state.p_history[i] = p_his;
    end;

    # set connection status to true
    xylem.state.connected = true;

    return nothing
);
