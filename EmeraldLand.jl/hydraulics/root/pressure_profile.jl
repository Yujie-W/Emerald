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
