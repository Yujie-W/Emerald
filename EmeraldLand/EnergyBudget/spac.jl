# This file contains function to run plant energy budgets due to water flow

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-30: add function to calculate the energy flow of the plant (due to water flow)
#
#######################################################################################################################################################################################################
"""

    spac_energy_flow!(spac::MultiLayerSPAC{FT}) where {FT}

Calculate the energy flows of the plant, given
- `spac` `MultiLayerSPAC` type SPAC

"""
function spac_energy_flow!(spac::MultiLayerSPAC{FT}) where {FT}
    root_energy_flows!(spac);
    junction_energy_flows!(spac);
    stem_energy_flows!(spac);
    leaf_energy_flows!(spac);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Oct-02: add function to run spac energy budget
#
#######################################################################################################################################################################################################
"""

    spac_energy_budget!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Calculate the energy budgets of the spac, given
- `config` `SPACConfiguration` type configuration
- `spac` `MultiLayerSPAC` type SPAC
- `δt` time step

"""
function spac_energy_budget!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}
    (; DIM_LAYER, DIM_ROOT) = config;
    (; AIR, BRANCHES, LEAVES, LEAVES_INDEX, ROOTS, TRUNK) = spac;

    # update the temperature for roots
    for i in 1:DIM_ROOT
        ROOTS[i].energy.state.Σe += ROOTS[i].energy.auxil.∂e∂t * δt;
    end;

    # update the temperature for trunk
    TRUNK.energy.state.Σe += TRUNK.energy.auxil.∂e∂t * δt;

    # update the temperature for branches and leaves
    for i in 1:DIM_LAYER
        BRANCHES[i].energy.state.Σe += BRANCHES[i].energy.auxil.∂e∂t * δt;
        LEAVES[i].energy.state.Σe += LEAVES[i].energy.auxil.∂e∂t * δt;
    end;

    return nothing
end;
