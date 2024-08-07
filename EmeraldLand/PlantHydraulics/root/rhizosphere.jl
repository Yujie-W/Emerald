# This file contains function to update rhizosphere pressure profile

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-23: add function rhizosphere_pressure_profile!
#
#######################################################################################################################################################################################################
"""

    rhizosphere_pressure_profile!(root::Root{FT}, soil::SoilLayer{FT}) where {FT}

Update the rhizosphere pressure profile, given
- `root` `Root` type struct
- `soil` `SoilLayer` type struct

"""
function rhizosphere_pressure_profile!(root::Root{FT}, soil::SoilLayer{FT}) where {FT}
    # compute the pressure at the end; of rhizosphere
    k_rhizo_max = root.rhizosphere.state.k_max * root.xylem.trait.area;
    f_st_soil = relative_surface_tension(soil.s_aux.t);
    f_vis_soil = relative_viscosity(soil.s_aux.t);
    f = flow_in(root);
    p = soil.s_aux.ψ;
    for _ in 1:5
        k = relative_soil_k(soil.trait.vc, true, p / f_st_soil) * k_rhizo_max * 5 / f_vis_soil;
        p -= f / k;
    end;
    root.rhizosphere.auxil.p_rhizo = p;

    return nothing
end;
