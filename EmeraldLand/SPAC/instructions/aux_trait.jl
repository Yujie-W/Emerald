# this file contains functions (methods) to update the trait-dependent auxiliary variables of the entire SPAC system

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-25: add t_aux! method for CanopyStructureTrait-dependent variables
#     2024-Feb-25: add t_aux! method for the combined MultiLayerCanopy
#     2024-Feb-27: add t_aux! method for the bulk SPAC system
#     2024-Feb-28: set x_bnds to 0 if LAI <= 0
#
#######################################################################################################################################################################################################
"""

    t_aux!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update the trait-dependent auxiliary variables for the SPAC system, given
- `config` SPAC configuration
- `spac` SPAC

"""
function t_aux! end;

t_aux!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    # update the t_aux for each of the field in the bulk spac system (the order should not matter)
    # the soil auxiliary variables
    for soil in spac.soils
        t_aux!(soil);
    end;

    # the air auxiliary variables
    for air in spac.airs
        t_aux!(air);
    end;

    # the canopy auxiliary variables
    t_aux!(config, spac.canopy);

    return nothing
);

t_aux!(soil::SoilLayer{FT}) where {FT} = (
    soil.t_aux.z = abs(soil.trait.zs[1] + soil.trait.zs[2]) / 2;
    soil.t_aux.δz = soil.trait.zs[1] - soil.trait.zs[2];

    return nothing
);

t_aux!(config::SPACConfiguration{FT}, canopy::MultiLayerCanopy{FT}) where {FT} = t_aux!(config, canopy.structure);

t_aux!(config::SPACConfiguration{FT}, canstr::CanopyStructure{FT}) where {FT} = (
    if canstr.trait.lai <= 0 && canstr.trait.sai <= 0
        canstr.t_aux.x_bnds .= 0;
    else
        canstr.t_aux.x_bnds .= [0; [sum(canstr.trait.δlai[1:i]) + sum(canstr.trait.δsai[1:i]) for i in eachindex(canstr.trait.δlai)]] ./ -(canstr.trait.lai + canstr.trait.sai);
    end;
    canopy_structure_aux!(config, canstr.trait, canstr.t_aux);

    return nothing
);

t_aux!(air::AirLayer{FT}) where {FT} = (
    air.t_aux.δz = air.trait.zs[2] - air.trait.zs[1];
    air.t_aux.z = (air.trait.zs[1] + air.trait.zs[2]) / 2;

    return nothing
);
