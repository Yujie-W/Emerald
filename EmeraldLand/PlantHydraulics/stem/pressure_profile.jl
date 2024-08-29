# This file contains functions related to stem pressure profile

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-25: add function stem_pressure_profile!
#
#######################################################################################################################################################################################################
"""

    stem_pressure_profile!(stem::Stem{FT}, p_dos::FT) where {FT}

Update the stem xylem pressure profile, given
- `stem` `Stem` type struct
- `p_dos` Downstream xylem pressure

"""
function stem_pressure_profile!(stem::Stem{FT}, p_dos::FT) where {FT}
    stem.xylem.auxil.pressure[1] = p_dos;
    xylem_pressure_profile!(stem.xylem, stem.energy.s_aux.t);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function stem_pressure_profiles!
#     2024-Feb-28: add LAI <= 0 control
#     2024-Aug-06: remove LAI control so that the function is always running to update the pressure for leaves to regrow
#
#######################################################################################################################################################################################################
"""

    stem_pressure_profiles!(spac::BulkSPAC{FT}) where {FT}

Set up stem pressure profile for trunk and branches, given
- `spac` `BulkSPAC` type struct

"""
function stem_pressure_profiles!(spac::BulkSPAC{FT}) where {FT}
    # run the pressure profile calculation for each stem layer only if LAI > 0
    branches = spac.plant.branches;
    junction = spac.plant.junction;
    trunk = spac.plant.trunk;

    stem_pressure_profile!(trunk, junction.s_aux.pressure);
    for stem in branches
        stem_pressure_profile!(stem, (trunk).xylem.auxil.pressure[end]);
    end;

    return nothing
end;
