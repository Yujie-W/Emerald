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
- Leaf (ind=NA for shaded, ind>1 for sunlit leaf)

"""
function empirical_equation end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for BallBerrySM using Leaf for shaded leaf
#     2022-Jul-07: add method for GentineSM using Leaf for shaded leaf
#     2022-Jul-07: add method for LeuningSM using Leaf for shaded leaf
#     2022-Jul-07: add method for MedlynSM using Leaf for shaded leaf
#     2022-Oct-20: add a max controller to make sure vpd is at least 1 Pa
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::BallBerrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT}
    empirical_equation(sm::GentineSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT}
    empirical_equation(sm::LeuningSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT}
    empirical_equation(sm::MedlynSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT}

Return the stomatal conductance computed from empirical model formulation for the shaded leaf of `Leaf`, given
- `sm` `BallBerrySM`, `GentineSM`, `LeuningSM`, or `MedlynSM` type model
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)

"""
empirical_equation(sm::BallBerrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;

    return G0 + β * G1 * air.auxil.ps[3] / saturation_vapor_pressure(air.auxil.t) * leaf.a_net_shaded * FT(1e-6) / leaf._p_CO₂_s_shaded * air.state.p_air
);

empirical_equation(sm::GentineSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;

    return G0 + β * G1 * leaf.a_net_shaded * FT(1e-6) / leaf._p_CO₂_i_shaded * air.state.p_air
);

empirical_equation(sm::LeuningSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; D0, G0, G1) = sm;

    _γ_s = (typeof(leaf.photosystem) <: C4VJP) ? 0 : leaf.photosystem.auxil.γ_star;
    _vpd = max(1, saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000) - air.auxil.ps[3]);

    return G0 + β * G1 / (1 + _vpd / D0) * leaf.a_net_shaded * FT(1e-6) / (leaf._p_CO₂_s_shaded - _γ_s) * air.state.p_air
);

empirical_equation(sm::MedlynSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;

    _vpd = max(1, saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000) - air.auxil.ps[3]);

    return G0 + FT(1.6) * (1 + β * G1 / sqrt(_vpd)) * leaf.a_net_shaded * FT(1e-6) / air.auxil.ps[2] * air.state.p_air
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method for BallBerrySM using Leaf for sunlit leaf
#     2022-Jul-07: add method for GentineSM using Leaf for sunlit leaf
#     2022-Jul-07: add method for LeuningSM using Leaf for sunlit leaf
#     2022-Jul-07: add method for MedlynSM using Leaf for sunlit leaf
#     2022-Oct-20: add a max controller to make sure vpd is at least 1 Pa
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm::BallBerrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT}
    empirical_equation(sm::GentineSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT}
    empirical_equation(sm::LeuningSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT}
    empirical_equation(sm::MedlynSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT}

Return the stomatal conductance computed from empirical model formulation for the sunlit leaf of `Leaf`, given
- `sm` `BallBerrySM`, `GentineSM`, `LeuningSM`, or `MedlynSM` type model
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `ind` Sunlit leaf index within the leaf angular distribution
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)

"""
empirical_equation(sm::BallBerrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;

    return G0 + β * G1 * air.auxil.ps[3] / saturation_vapor_pressure(air.auxil.t) * leaf.flux.auxil.a_n_sunlit[ind] * FT(1e-6) / leaf._p_CO₂_s_sunlit[ind] * air.state.p_air
);

empirical_equation(sm::GentineSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;

    return G0 + β * G1 * leaf.flux.auxil.a_n_sunlit[ind] * FT(1e-6) / leaf._p_CO₂_i_sunlit[ind] * air.state.p_air
);

empirical_equation(sm::LeuningSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT} = (
    (; D0, G0, G1) = sm;

    _γ_s = (typeof(leaf.photosystem) <: C4VJP) ? 0 : leaf.photosystem.auxil.γ_star;
    _vpd = max(1, saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000) - air.auxil.ps[3]);

    return G0 + β * G1 / (1 + _vpd / D0) * leaf.flux.auxil.a_n_sunlit[ind] * FT(1e-6) / (leaf._p_CO₂_s_sunlit[ind] - _γ_s) * air.state.p_air
);

empirical_equation(sm::MedlynSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;

    _vpd = max(1, saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000) - air.auxil.ps[3]);

    return G0 + FT(1.6) * (1 + β * G1 / sqrt(_vpd)) * leaf.flux.auxil.a_n_sunlit[ind] * FT(1e-6) / air.auxil.ps[2] * air.state.p_air
);
