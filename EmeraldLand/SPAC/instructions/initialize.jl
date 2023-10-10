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
    soil.state.ns[4] = air.state.p_air * 0.79 * soil.auxil.δz * δθ / rt;
    soil.state.ns[5] = air.state.p_air * 0.209 * soil.auxil.δz * δθ / rt;
    soil.auxil.cp = heat_capacitance(soil);
    soil.state.Σe = soil.auxil.cp * soil.auxil.t;

    return nothing
);

initialize_struct!(root::Root{FT}) where {FT} = (
    root.xylem.state.v_storage .= (root.xylem.state.v_max * root.xylem.state.area * root.xylem.state.l) / length(root.xylem.state.v_storage);
    root.energy.auxil.cp = sum(root.xylem.state.v_storage) * CP_L_MOL(FT) + (root.xylem.state.cp * root.xylem.state.area * root.xylem.state. l);
    root.energy.state.Σe = root.energy.auxil.cp * root.energy.auxil.t;

    return nothing
);

initialize_struct!(stem::Stem{FT}) where {FT} = (
    stem.xylem.state.v_storage .= (stem.xylem.state.v_max * stem.xylem.state.area * stem.xylem.state.l) / length(stem.xylem.state.v_storage);
    stem.energy.auxil.cp = sum(stem.xylem.state.v_storage) * CP_L_MOL(FT) + (stem.xylem.state.cp * stem.xylem.state.area * stem.xylem.state. l);
    stem.energy.state.Σe = stem.energy.auxil.cp * stem.energy.auxil.t;

    return nothing
);

initialize_struct!(leaf::Leaf{FT}) where {FT} = (
    leaf.xylem.state.cp = 1780;
    leaf.xylem.state.k_max = 0.04;
    leaf.capacitor.state.v_storage = leaf.capacitor.state.v_max * leaf.xylem.state.area;
    leaf.energy.auxil.cp = leaf.capacitor.state.v_storage * CP_L_MOL(FT) + leaf.bio.state.lma * 10 * leaf.xylem.state.area * leaf.xylem.state.cp;
    leaf.energy.state.Σe = leaf.energy.auxil.cp * leaf.energy.auxil.t;

    return nothing
);

initialize_struct!(air::AirLayer{FT}) where {FT} = (
    air.auxil.ps[2] = air.auxil.f_CO₂ * air.state.p_air * 1e-6;
    air.auxil.ps[4] = 0.79 * air.state.p_air;
    air.auxil.ps[5] = 0.209 * air.state.p_air;

    air.state.ns[2] = air.auxil.f_CO₂ * air.state.p_air * 1e-6 * air.auxil.δz / (GAS_R(FT) * air.auxil.t);
    air.state.ns[3] = air.auxil.ps[3] * air.auxil.δz / (GAS_R(FT) * air.auxil.t);
    air.state.ns[4] = 0.79 * air.state.p_air * air.auxil.δz / (GAS_R(FT) * air.auxil.t);
    air.state.ns[5] = 0.209 * air.state.p_air * air.auxil.δz / (GAS_R(FT) * air.auxil.t);

    air.state.Σe = heat_capacitance(air) * air.auxil.t;

    return nothing
);
