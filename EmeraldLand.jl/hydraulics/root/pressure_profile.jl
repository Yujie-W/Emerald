# This file contains functions related to root pressure profile

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-23: add function root_pressure_profile!
#
#######################################################################################################################################################################################################
"""

    root_pressure_profile!(root::Root{FT}, soil::SoilLayer{FT}) where {FT}

Update the rhizosphere and root xylem pressure profile, given
- `root` `Root` type struct
- `soil` `SoilLayer` type struct

"""
function root_pressure_profile!(root::Root{FT}, soil::SoilLayer{FT}) where {FT}
    rhizosphere_pressure_profile!(root, soil);
    root.xylem.auxil.pressure[1] = root.rhizosphere.auxil.p_rhizo;
    xylem_pressure_profile!(root.xylem, root.energy.auxil.t);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function root_pressure_profiles!
#
#######################################################################################################################################################################################################
"""

    root_pressure_profiles!(spac::MultiLayerSPAC{FT}) where {FT}

Set up root pressure profile for each root, given
- `spac` `MultiLayerSPAC` type struct

"""
function root_pressure_profiles!(spac::MultiLayerSPAC{FT}) where {FT}
    (; ROOTS, ROOTS_INDEX, SOIL) = spac;

    for i in eachindex(ROOTS)
        root_pressure_profile!(ROOTS[i].NS, SOIL.LAYERS[ROOTS_INDEX[i]]);
    end;

    return nothing
end;
