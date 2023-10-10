# This file contains function to update leaf water pressure profile along the capacitor

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-26: add function to update the extraxylary pressure profile
#
#######################################################################################################################################################################################################
"""

    extraxylary_pressure_profile!(leaf::Leaf{FT}) where {FT}

Update the extraxylary pressure profile, given
- `leaf` `Leaf` type struct

"""
function extraxylary_pressure_profile! end;

extraxylary_pressure_profile!(leaf::Leaf{FT}) where {FT} = extraxylary_pressure_profile!(leaf.xylem.state, leaf.xylem.auxil, leaf.capacitor.state, leaf.capacitor.auxil, leaf.energy.auxil.t);

extraxylary_pressure_profile!(
            x_state::XylemHydraulicsState{FT},
            x_auxil::XylemHydraulicsAuxilNSS{FT},
            c_state::ExtraXylemCapacitorState{FT},
            c_auxil::ExtraXylemCapacitorAuxil{FT},
            t::FT) where {FT} = (
    flow = flow_out(x_auxil) + c_auxil.flow;
    c_auxil.p_leaf = x_auxil.pressure[end;] - flow / c_auxil.k;

    return nothing
);

extraxylary_pressure_profile!(
            x_state::XylemHydraulicsState{FT},
            x_auxil::XylemHydraulicsAuxilSS{FT},
            c_state::ExtraXylemCapacitorState{FT},
            c_auxil::ExtraXylemCapacitorAuxil{FT},
            t::FT) where {FT} = (
    k_max = x_state.area * c_state.k_max;
    f_st = relative_surface_tension(t);
    f_vis = relative_viscosity(t);
    flow = flow_out(x_auxil);
    k = relative_xylem_k(c_state.vc, x_auxil.pressure[end;] / f_st) / f_vis * k_max;
    c_auxil.p_leaf = x_auxil.pressure[end;] - flow / k;

    return nothing
);
