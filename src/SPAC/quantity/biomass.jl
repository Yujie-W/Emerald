# This file is meant to save the biomass related parameters

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Sep-09: add function to return the total sap volume of the plant (to compute mortality rate)
#
#######################################################################################################################################################################################################
"""

    SAP_VOLUME(spac::BulkSPAC{FT}) where {FT}

Return the total sap volume of the plant, given
- `spac` the SPAC struct

"""
function SAP_VOLUME(spac::BulkSPAC{FT}) where {FT}
    sap_vol::FT = 0;

    # loop through all the roots
    for r in spac.plant.roots
        sap_vol += r.xylem.state.asap * r.xylem.trait.l;
    end;

    # the trunk
    sap_vol += spac.plant.trunk.xylem.state.asap * spac.plant.trunk.xylem.trait.l;

    # loop through all the branches
    for s in spac.plant.branches
        sap_vol += s.xylem.state.asap * s.xylem.trait.l;
    end;

    return sap_vol
end;
