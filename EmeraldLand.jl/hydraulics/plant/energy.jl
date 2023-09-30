# This file contains function to run plant energy budgets due to water flow

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-30: add function to calculate the energy flow of the plant (due to water flow)
#
#######################################################################################################################################################################################################
"""

    plant_energy_flow!(spac::MultiLayerSPAC{FT}) where {FT}

Calculate the energy flows of the plant, given
- `spac` `MultiLayerSPAC` type SPAC

"""
function plant_energy_flow!(spac::MultiLayerSPAC{FT}) where {FT}
    root_energy_flows!(spac);
    junction_energy_flows!(spac);
    stem_energy_flows!(spac);
    leaf_energy_flows!(spac);

    return nothing
end;
