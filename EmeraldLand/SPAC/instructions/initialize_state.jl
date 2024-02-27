# This file contains function to initialize the structs used in the SPAC model

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Oct-09: add function to initialize SPAC states
#     2024-Feb-23: rename function to initialize_states! to be more consistent with its aim
#     2024-Feb-27: rename function to initialize_energy_states! to be more consistent with its aim
#
#######################################################################################################################################################################################################
"""

Initialize the energy related state variables of the SPAC structs.

"""
function initialize_energy_states! end;

initialize_energy_states!(soil::SoilLayer{FT}, air::AirLayer{FT}) where {FT} = (
    δθ = max(0, soil.trait.vc.Θ_SAT - soil.state.θ);
    rt = GAS_R(FT) * soil.s_aux.t;
    soil.state.ns[3] = saturation_vapor_pressure(soil.s_aux.t) * soil.t_aux.δz * (δθ + FT(0.01)) / rt;
    soil.state.ns[4] = air.state.p_air * F_N₂(FT) * soil.t_aux.δz * δθ / rt;
    soil.state.ns[5] = air.state.p_air * F_O₂(FT) * soil.t_aux.δz * δθ / rt;
    soil.s_aux.cp = heat_capacitance(soil);
    soil.state.Σe = soil.s_aux.cp * soil.s_aux.t;

    return nothing
);

initialize_energy_states!(root::Root{FT}) where {FT} = (
    root.xylem.state.v_storage .= (root.xylem.trait.v_max * root.xylem.trait.area * root.xylem.trait.l) / length(root.xylem.state.v_storage);
    root.energy.s_aux.cp = heat_capacitance(root);
    root.energy.state.Σe = root.energy.s_aux.cp * root.energy.s_aux.t;

    return nothing
);

initialize_energy_states!(stem::Stem{FT}) where {FT} = (
    stem.xylem.state.v_storage .= (stem.xylem.trait.v_max * stem.xylem.trait.area * stem.xylem.trait.l) / length(stem.xylem.state.v_storage);
    stem.energy.s_aux.cp = heat_capacitance(stem);
    stem.energy.state.Σe = stem.energy.s_aux.cp * stem.energy.s_aux.t;

    return nothing
);

initialize_energy_states!(leaf::Leaf{FT}) where {FT} = (
    leaf.capacitor.state.v_storage = leaf.capacitor.trait.v_max;
    leaf.energy.s_aux.cp = heat_capacitance(leaf);
    leaf.energy.state.Σe = leaf.energy.s_aux.cp * leaf.energy.s_aux.t;

    return nothing
);

initialize_energy_states!(air::AirLayer{FT}) where {FT} = (
    air.s_aux.ps[2] = air.s_aux.f_CO₂ * air.state.p_air * 1e-6;
    air.s_aux.ps[4] = F_N₂(FT) * air.state.p_air;
    air.s_aux.ps[5] = F_O₂(FT) * air.state.p_air;

    for i in 1:5
        air.state.ns[i] = air.s_aux.ps[i] * air.t_aux.δz / (GAS_R(FT) * air.s_aux.t);
    end;

    air.state.Σe = heat_capacitance(air) * air.s_aux.t;

    return nothing
);

initialize_energy_states!(spac::BulkSPAC{FT}) where {FT} = (
    airs = spac.airs;
    branches = spac.plant.branches;
    leaves = spac.plant.leaves;
    roots = spac.plant.roots;
    soils = spac.soils;
    trunk = spac.plant.trunk;

    # make sure soil energy is correctly scaled with temperature and soil water content
    for soil in soils
        initialize_energy_states!(soil, airs[1]);
    end;

    # make sure the root energy is correctly scaled with temperature
    for root in roots
        initialize_energy_states!(root);
    end;

    # make sure the stem energy is correctly scaled with temperature
    initialize_energy_states!(trunk);
    for stem in branches
        initialize_energy_states!(stem);
    end;

    # make sure leaf area index setup and energy are correct
    for leaf in leaves
        initialize_energy_states!(leaf);
    end;

    # make sure air layers are correctly initialized
    for air in airs
        initialize_energy_states!(air);
    end;

    return nothing
);
