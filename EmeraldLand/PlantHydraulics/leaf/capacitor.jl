# This file contains function to update leaf water pressure profile along the capacitor

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-26: add function to update the extraxylary pressure profile
#     2024-Feb-28: p_leaf now is a state (see the comment in the state struct for the reason)
#
#######################################################################################################################################################################################################
"""

    extraxylary_pressure_profile!(leaf::Leaf{FT}) where {FT}

Update the extraxylary pressure profile, given
- `leaf` `Leaf` type struct

"""
function extraxylary_pressure_profile! end;

extraxylary_pressure_profile!(leaf::Leaf{FT}) where {FT} =
    extraxylary_pressure_profile!(leaf.xylem.trait, leaf.xylem.auxil, leaf.capacitor.trait, leaf.capacitor.state, leaf.capacitor.auxil, leaf.energy.s_aux.t);

extraxylary_pressure_profile!(
            x_trait::XylemHydraulicsTrait{FT},
            x_auxil::XylemHydraulicsAuxilNSS{FT},
            c_trait::ExtraXylemCapacitorTrait{FT},
            c_state::ExtraXylemCapacitorState{FT},
            c_auxil::ExtraXylemCapacitorAuxil{FT},
            t::FT) where {FT} = (
    flow = flow_out(x_auxil) + c_auxil.flow;
    c_state.p_leaf = x_auxil.pressure[end] - flow / c_auxil.k;

    return nothing
);

extraxylary_pressure_profile!(
            x_trait::XylemHydraulicsTrait{FT},
            x_auxil::XylemHydraulicsAuxilSS{FT},
            c_trait::ExtraXylemCapacitorTrait{FT},
            c_state::ExtraXylemCapacitorState{FT},
            c_auxil::ExtraXylemCapacitorAuxil{FT},
            t::FT) where {FT} = (
    k_max = x_trait.area * c_trait.k_max;
    f_st = relative_surface_tension(t);
    f_vis = relative_viscosity(t);
    flow = flow_out(x_auxil);
    k = relative_xylem_k(c_trait.vc, x_auxil.pressure[end] / f_st) / f_vis * k_max;
    c_state.p_leaf = x_auxil.pressure[end] - flow / k;

    return nothing
);
