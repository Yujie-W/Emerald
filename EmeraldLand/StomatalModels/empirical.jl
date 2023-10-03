#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-01: migrate function from older version
#     2022-Jul-01: rename the function from stomatal_conductance to empirical_equation
#
#######################################################################################################################################################################################################
"""
This function returns the stomatal conductance computed from empirical stomatal models. This is not the solution! Supported methods are for
- Leaf
- Leaves2D (ind=NA for shaded, ind>1 for sunlit leaves)

"""
function empirical_equation end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for BallBerrySM using Leaves2D for shaded leaves
#     2022-Jul-07: add method for GentineSM using Leaves2D for shaded leaves
#     2022-Jul-07: add method for LeuningSM using Leaves2D for shaded leaves
#     2022-Jul-07: add method for MedlynSM using Leaves2D for shaded leaves
#     2022-Oct-20: add a max controller to make sure vpd is at least 1 Pa
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT}
    empirical_equation(sm::GentineSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT}
    empirical_equation(sm::LeuningSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT}
    empirical_equation(sm::MedlynSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT}

Return the stomatal conductance computed from empirical model formulation for the shaded leaves of `Leaves2D`, given
- `sm` `BallBerrySM`, `GentineSM`, `LeuningSM`, or `MedlynSM` type model
- `leaves` `Leaves2D` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)

"""
empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;
    (; P_AIR) = air;

    return G0 + β * G1 * air.p_H₂O / saturation_vapor_pressure(air.t) * leaves.a_net_shaded * FT(1e-6) / leaves._p_CO₂_s_shaded * P_AIR
);

empirical_equation(sm::GentineSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;
    (; P_AIR) = air;

    return G0 + β * G1 * leaves.a_net_shaded * FT(1e-6) / leaves._p_CO₂_i_shaded * P_AIR
);

empirical_equation(sm::LeuningSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; D0, G0, G1) = sm;
    (; P_AIR) = air;

    _γ_s = (typeof(leaves.NS.photosystem) <: C4VJP) ? 0 : leaves.PSM._γ_star;
    _vpd = max(1, saturation_vapor_pressure(leaves.NS.energy.auxil.t, leaves.NS.capacitor.auxil.p_leaf * 1000000) - air.p_H₂O);

    return G0 + β * G1 / (1 + _vpd / D0) * leaves.a_net_shaded * FT(1e-6) / (leaves._p_CO₂_s_shaded - _γ_s) * P_AIR
);

empirical_equation(sm::MedlynSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;
    (; P_AIR) = air;

    _vpd = max(1, saturation_vapor_pressure(leaves.NS.energy.auxil.t, leaves.NS.capacitor.auxil.p_leaf * 1000000) - air.p_H₂O);

    return G0 + FT(1.6) * (1 + β * G1 / sqrt(_vpd)) * leaves.a_net_shaded * FT(1e-6) / air.p_CO₂ * P_AIR
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for BallBerrySM using Leaves2D for sunlit leaves
#     2022-Jul-07: add method for GentineSM using Leaves2D for sunlit leaves
#     2022-Jul-07: add method for LeuningSM using Leaves2D for sunlit leaves
#     2022-Jul-07: add method for MedlynSM using Leaves2D for sunlit leaves
#     2022-Oct-20: add a max controller to make sure vpd is at least 1 Pa
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT}
    empirical_equation(sm::GentineSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT}
    empirical_equation(sm::LeuningSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT}
    empirical_equation(sm::MedlynSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT}

Return the stomatal conductance computed from empirical model formulation for the sunlit leaves of `Leaves2D`, given
- `sm` `BallBerrySM`, `GentineSM`, `LeuningSM`, or `MedlynSM` type model
- `leaves` `Leaves2D` type struct
- `air` `AirLayer` type environmental conditions
- `ind` Sunlit leaf index within the leaf angular distribution
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)

"""
empirical_equation(sm::BallBerrySM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;
    (; P_AIR) = air;

    return G0 + β * G1 * air.p_H₂O / saturation_vapor_pressure(air.t) * leaves.a_net_sunlit[ind] * FT(1e-6) / leaves._p_CO₂_s_sunlit[ind] * P_AIR
);

empirical_equation(sm::GentineSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;
    (; P_AIR) = air;

    return G0 + β * G1 * leaves.a_net_sunlit[ind] * FT(1e-6) / leaves._p_CO₂_i_sunlit[ind] * P_AIR
);

empirical_equation(sm::LeuningSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT} = (
    (; D0, G0, G1) = sm;
    (; P_AIR) = air;

    _γ_s = (typeof(leaves.NS.photosystem) <: C4VJP) ? 0 : leaves.PSM._γ_star;
    _vpd = max(1, saturation_vapor_pressure(leaves.NS.energy.auxil.t, leaves.NS.capacitor.auxil.p_leaf * 1000000) - air.p_H₂O);

    return G0 + β * G1 / (1 + _vpd / D0) * leaves.a_net_sunlit[ind] * FT(1e-6) / (leaves._p_CO₂_s_sunlit[ind] - _γ_s) * P_AIR
);

empirical_equation(sm::MedlynSM{FT}, leaves::Leaves2D{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;
    (; P_AIR) = air;

    _vpd = max(1, saturation_vapor_pressure(leaves.NS.energy.auxil.t, leaves.NS.capacitor.auxil.p_leaf * 1000000) - air.p_H₂O);

    return G0 + FT(1.6) * (1 + β * G1 / sqrt(_vpd)) * leaves.a_net_sunlit[ind] * FT(1e-6) / air.p_CO₂ * P_AIR
);
