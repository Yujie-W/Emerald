# This file contains function to run plant energy budgets due to water flow

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-30: add function to calculate the energy flow of the plant (due to water flow)
#     2023-Oct-07: add soil energy flow
#
#######################################################################################################################################################################################################
"""

    spac_energy_flow!(spac::MultiLayerSPAC{FT}) where {FT}

Calculate the energy flows of the plant, given
- `spac` `MultiLayerSPAC` type SPAC

"""
function spac_energy_flow!(spac::MultiLayerSPAC{FT}) where {FT}
    soil_energy_flow!(spac);
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
#     2023-Oct-07: add soil energy budget
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
    (; BRANCHES, LEAVES, ROOTS, SOILS, TRUNK) = spac;

    # update the temperature for soil
    for soil in SOILS
        # water mass and energy flow
        soil.state.Σe += soil.auxil.∂e∂t * δt / soil.auxil.δz;

        # soil water condensation or evaporation
        soil.state.Σe += soil.auxil.n_con * M_H₂O(FT) * latent_heat_vapor(soil.auxil.t) / soil.auxil.δz;
    end;

    # update the energy loss related to surface runoff
    if SOILS[1].state.θ > SOILS[1].state.vc.Θ_SAT
        cp = heat_capacitance(SOILS[1]; runoff = SOILS[1].auxil.runoff);
        t  = SOILS[1].state.Σe / cp;
        SOILS[1].state.Σe -= SOILS[1].auxil.runoff / SOILS[1].auxil.δz * CP_L_MOL(FT) * t;
    end;

    # update the temperature for roots
    for root in ROOTS
        root.energy.state.Σe += root.energy.auxil.∂e∂t * δt;
    end;

    # update the temperature for trunk
    TRUNK.energy.state.Σe += TRUNK.energy.auxil.∂e∂t * δt;

    # update the temperature for branches
    for stem in BRANCHES
        stem.energy.state.Σe += stem.energy.auxil.∂e∂t * δt;
    end;

    # update the temperature for leaves
    for leaf in LEAVES
        leaf.energy.state.Σe += leaf.energy.auxil.∂e∂t * δt;
    end;

    return nothing
end;
