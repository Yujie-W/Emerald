# This file contains functions related to root pressure profile

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-23: add function root_pressure_profile!
#     2023-Sep-30: if root is diconnected, update pressure from junction to root xylem (excluding rhizosphere)
#
#######################################################################################################################################################################################################
"""

    root_pressure_profile!(soil::SoilLayer{FT}, root::Root{FT}, junction::JunctionCapacitor{FT}) where {FT}

Update the rhizosphere and root xylem pressure profile, given
- `soil` `SoilLayer` type struct
- `root` `Root` type struct
- `junction` `JunctionCapacitor` type struct

"""
function root_pressure_profile!(soil::SoilLayer{FT}, root::Root{FT}, junction::JunctionCapacitor{FT}) where {FT}
    # if root is connected, update pressure from soil to rhizosphere and root xylem
    # else, update pressure from junction to root xylem (excluding rhizosphere)
    if root.xylem.auxil.connected
        rhizosphere_pressure_profile!(root, soil);
        root.xylem.auxil.pressure[1] = root.rhizosphere.auxil.p_rhizo;
        xylem_pressure_profile!(root.xylem, root.energy.auxil.t);
    else
        root.xylem.auxil.pressure[end] = junction.auxil.pressure;
        xylem_pressure_profile!(root.xylem, root.energy.auxil.t, true);
    end;

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

    root_pressure_profiles!(spac::BulkSPAC{FT}) where {FT}

Set up root pressure profile for each root, given
- `spac` `BulkSPAC` type struct

"""
function root_pressure_profiles!(spac::BulkSPAC{FT}) where {FT}
    (; JUNCTION, ROOTS, ROOTS_INDEX, SOILS) = spac;

    for i in eachindex(ROOTS)
        root_pressure_profile!(SOILS[ROOTS_INDEX[i]], ROOTS[i], JUNCTION);
    end;

    return nothing
end;
