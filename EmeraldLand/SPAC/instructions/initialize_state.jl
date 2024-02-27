# This file contains function to initialize the structs used in the SPAC model

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Oct-09: add function to initialize SPAC states
#     2024-Feb-23: rename function to initialize_states! to be more consistent with its aim
#
#######################################################################################################################################################################################################
"""

Initialize the energy and state variables of the SPAC structs.

"""
function initialize_states! end;

initialize_states!(soil::SoilLayer{FT}, air::AirLayer{FT}) where {FT} = (
    δθ = max(0, soil.trait.vc.Θ_SAT - soil.state.θ);
    rt = GAS_R(FT) * soil.s_aux.t;
    soil.state.ns[3] = saturation_vapor_pressure(soil.s_aux.t, soil.s_aux.ψ * 1000000) * soil.t_aux.δz * (δθ + FT(0.01)) / rt;
    soil.state.ns[4] = air.state.p_air * F_N₂(FT) * soil.t_aux.δz * δθ / rt;
    soil.state.ns[5] = air.state.p_air * F_O₂(FT) * soil.t_aux.δz * δθ / rt;
    soil.s_aux.cp = heat_capacitance(soil);
    soil.state.Σe = soil.s_aux.cp * soil.s_aux.t;

    return nothing
);

initialize_states!(root::Root{FT}) where {FT} = (
    root.xylem.state.v_storage .= (root.xylem.trait.v_max * root.xylem.trait.area * root.xylem.trait.l) / length(root.xylem.state.v_storage);
    root.energy.s_aux.cp = heat_capacitance(root);
    root.energy.state.Σe = root.energy.s_aux.cp * root.energy.s_aux.t;

    return nothing
);

initialize_states!(stem::Stem{FT}) where {FT} = (
    stem.xylem.state.v_storage .= (stem.xylem.trait.v_max * stem.xylem.trait.area * stem.xylem.trait.l) / length(stem.xylem.state.v_storage);
    stem.energy.s_aux.cp = heat_capacitance(stem);
    stem.energy.state.Σe = stem.energy.s_aux.cp * stem.energy.s_aux.t;

    return nothing
);

initialize_states!(leaf::Leaf{FT}; k_sla::Number = 0.04) where {FT} = (
    leaf.xylem.trait.cp = 1780;
    leaf.xylem.trait.k_max = k_sla;
    leaf.capacitor.state.v_storage = leaf.capacitor.trait.v_max;
    leaf.energy.s_aux.cp = heat_capacitance(leaf);
    leaf.energy.state.Σe = leaf.energy.s_aux.cp * leaf.energy.s_aux.t;

    return nothing
);

initialize_states!(air::AirLayer{FT}) where {FT} = (
    air.s_aux.ps[2] = air.s_aux.f_CO₂ * air.state.p_air * 1e-6;
    air.s_aux.ps[4] = F_N₂(FT) * air.state.p_air;
    air.s_aux.ps[5] = F_O₂(FT) * air.state.p_air;

    for i in 1:5
        air.state.ns[i] = air.s_aux.ps[i] * air.t_aux.δz / (GAS_R(FT) * air.s_aux.t);
    end;

    air.state.Σe = heat_capacitance(air) * air.s_aux.t;

    return nothing
);

initialize_states!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    airs = spac.airs;
    branches = spac.plant.branches;
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    roots = spac.plant.roots;
    soils = spac.soils;
    sbulk = spac.soil_bulk;
    trunk = spac.plant.trunk;
    n_layer = length(leaves);

    # make sure soil energy is correctly scaled with temperature and soil water content
    for soil in soils
        initialize_states!(soil, airs[1]);
    end;

    # make sure the root energy is correctly scaled with temperature
    for root in roots
        initialize_states!(root);
    end;

    # make sure the stem energy is correctly scaled with temperature
    initialize_states!(trunk);
    for stem in branches
        initialize_states!(stem);
    end;

    # make sure leaf area index setup and energy are correct
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaves[ilf].xylem.trait.area = sbulk.trait.area * can_str.trait.δlai[irt];
        initialize_states!(leaves[ilf]);
    end;

    # make sure air layers are correctly initialized
    for air in airs
        initialize_states!(air);
    end;

    # initialize stomatal conductance
    stomatal_conductance!(spac, FT(0));

    return nothing
);
