# this file contains functions (methods) to update the trait-dependent auxilary variables of the entire SPAC system

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-27: add t_aux! method for the bulk SPAC system
#
#######################################################################################################################################################################################################
t_aux!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    # update the t_aux for each of the field in the bulk spac system (the order should not matter)
    # the soil auxilary variables
    for soil in spac.soils
        t_aux!(soil);
    end;

    # the air auxilary variables
    for air in spac.airs
        t_aux!(air);
    end;

    # the canopy auxilary variables
    t_aux!(config, spac.canopy);

    return nothing
);
