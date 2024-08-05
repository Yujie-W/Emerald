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

    spac_energy_flow!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Calculate the energy flows of the plant, given
- `config` `SPACConfiguration` type configuration
- `spac` `BulkSPAC` type SPAC

"""
function spac_energy_flow!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    soil_energy_flow!(spac);
    root_energy_flows!(spac);
    junction_energy_flows!(spac);
    stem_energy_flows!(spac);
    leaf_energy_flows!(config, spac);

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

    spac_energy_budget!(spac::BulkSPAC{FT}, δt::FT) where {FT}

Calculate the energy budgets of the spac, given
- `spac` `BulkSPAC` type SPAC
- `δt` time step

"""
function spac_energy_budget!(spac::BulkSPAC{FT}, δt::FT) where {FT}
    branches = spac.plant.branches;
    junction = spac.plant.junction;
    leaves = spac.plant.leaves;
    roots = spac.plant.roots;
    soils = spac.soils;
    trunk = spac.plant.trunk;

    # update the temperature for soil
    for soil in soils
        # water mass and energy flow
        soil.state.Σe += soil.auxil.∂e∂t * δt / soil.t_aux.δz;

        # soil water condensation or evaporation
        soil.state.Σe += soil.auxil.n_con * M_H₂O(FT) * latent_heat_vapor(soil.s_aux.t) / soil.t_aux.δz;
    end;

    # update the energy loss related to surface runoff
    top_soil = soils[1];
    if top_soil.state.θ > top_soil.trait.vc.Θ_SAT
        cp = heat_capacitance(top_soil; runoff = top_soil.auxil.runoff);
        t  = top_soil.state.Σe / cp;
        top_soil.state.Σe -= top_soil.auxil.runoff / top_soil.t_aux.δz * CP_L_MOL(FT) * t;
    end;

    # update the temperature for roots
    for root in roots
        root.energy.state.Σe += root.energy.auxil.∂e∂t * δt;
    end;

    # update the temperature for junction
    junction.state.Σe += junction.auxil.∂e∂t * δt;

    # update the temperature for trunk
    trunk.energy.state.Σe += trunk.energy.auxil.∂e∂t * δt;

    # update the temperature for branches
    for stem in branches
        stem.energy.state.Σe += stem.energy.auxil.∂e∂t * δt;
    end;

    # update the temperature for leaves
    for leaf in leaves
        leaf.energy.state.Σe += leaf.energy.auxil.∂e∂t * δt;
    end;

    return nothing
end;
