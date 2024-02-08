# This file contains function to initialize the structs used in the SPAC model

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Oct-09: add function to initialize SPAC structs
#
#######################################################################################################################################################################################################
"""

Initialize the energy and state variables of the SPAC structs.

"""
function initialize_struct! end;

initialize_struct!(soil::SoilLayer{FT}, air::AirLayer{FT}) where {FT} = (
    δθ = max(0, soil.state.vc.Θ_SAT - soil.state.θ);
    rt = GAS_R(FT) * soil.auxil.t;
    soil.state.ns[3] = saturation_vapor_pressure(soil.auxil.t, soil.auxil.ψ * 1000000) * soil.auxil.δz * (δθ + FT(0.01)) / rt;
    soil.state.ns[4] = air.state.p_air * F_N₂(FT) * soil.auxil.δz * δθ / rt;
    soil.state.ns[5] = air.state.p_air * F_O₂(FT) * soil.auxil.δz * δθ / rt;
    soil.auxil.cp = heat_capacitance(soil);
    soil.state.Σe = soil.auxil.cp * soil.auxil.t;

    return nothing
);

initialize_struct!(root::Root{FT}) where {FT} = (
    root.xylem.state.v_storage .= (root.xylem.state.v_max * root.xylem.state.area * root.xylem.state.l) / length(root.xylem.state.v_storage);
    root.energy.auxil.cp = heat_capacitance(root);
    root.energy.state.Σe = root.energy.auxil.cp * root.energy.auxil.t;

    return nothing
);

initialize_struct!(stem::Stem{FT}) where {FT} = (
    stem.xylem.state.v_storage .= (stem.xylem.state.v_max * stem.xylem.state.area * stem.xylem.state.l) / length(stem.xylem.state.v_storage);
    stem.energy.auxil.cp = heat_capacitance(stem);
    stem.energy.state.Σe = stem.energy.auxil.cp * stem.energy.auxil.t;

    return nothing
);

initialize_struct!(leaf::Leaf{FT}; k_sla::Number = 0.04) where {FT} = (
    leaf.xylem.state.cp = 1780;
    leaf.xylem.state.k_max = k_sla;
    leaf.capacitor.state.v_storage = leaf.capacitor.state.v_max;
    leaf.energy.auxil.cp = heat_capacitance(leaf);
    leaf.energy.state.Σe = leaf.energy.auxil.cp * leaf.energy.auxil.t;

    return nothing
);

initialize_struct!(air::AirLayer{FT}) where {FT} = (
    air.auxil.ps[2] = air.auxil.f_CO₂ * air.state.p_air * 1e-6;
    air.auxil.ps[4] = F_N₂(FT) * air.state.p_air;
    air.auxil.ps[5] = F_O₂(FT) * air.state.p_air;

    for i in 1:5
        air.state.ns[i] = air.auxil.ps[i] * air.auxil.δz / (GAS_R(FT) * air.auxil.t);
    end;

    air.state.Σe = heat_capacitance(air) * air.auxil.t;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-27: add function to initialize SPAC
#     2022-Jun-27: add leaf area controller to make sure soil and leaf areas are consistent with leaf area index
#     2023-Mar-27: initialize soil and leaf e as well (because T, SWC may be changed)
#     2023-Jun-12: initialize soil trace gas as well
#     2023-Jun-13: update N₂ and O₂ based on soil water content
#     2023-Jun-13: add soil gas energy into soil e
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Oct-07: add 0.01 to the water vapor volume per soil layer
#     2023-Oct-09: add root and stem initialization in the initialization of SPAC
#     2023-Oct-17: update step and subtep auxils during initialization
#     2023-Oct-18: initialize leaf inclination angles and canopy structure during initialization
#
#######################################################################################################################################################################################################
"""

    initialize!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Initialize the SPAC, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` SPAC

"""
function initialize! end;

initialize!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
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
        initialize_struct!(soil, airs[1]);
    end;

    # make sure the root energy is correctly scaled with temperature
    for root in roots
        initialize_struct!(root);
    end;

    # make sure the stem energy is correctly scaled with temperature
    initialize_struct!(trunk);
    for stem in branches
        initialize_struct!(stem);
    end;

    # make sure leaf area index setup and energy are correct
    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        leaves[ilf].xylem.state.area = sbulk.state.area * can_str.state.δlai[irt];
        initialize_struct!(leaves[ilf]);
    end;

    # make sure air layers are correctly initialized
    for air in airs
        initialize_struct!(air);
    end;

    # make sure tha auxilary variables are correctly initialized
    update_step_auxils!(spac);
    update_substep_auxils!(spac);

    # initialize leaf level spectra
    plant_leaf_spectra!(config, spac);

    # initialize canopy structure
    inclination_angles!(config, spac);
    canopy_structure!(config, spac);

    # initialize stomatal conductance
    stomatal_conductance!(spac, FT(0));

    return nothing
);
