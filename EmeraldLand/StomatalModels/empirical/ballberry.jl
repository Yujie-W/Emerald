# This file contains the function for Ball Berry stomatal model

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-01: migrate function from older version
#
#######################################################################################################################################################################################################
"""

    empirical_equation(sm, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT}
    empirical_equation(sm, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT}

Return the stomatal conductance computed from empirical model formulation for the shaded leaf of `Leaf`, given
- `sm` `BallBerrySM`, `GentineSM`, `LeuningSM`, or `MedlynSM` type empirical stomatal model
- `leaf` `Leaf` type struct
- `air` `AirLayer` type environmental conditions
- `β` Tuning factor for G1 (must be 1 if tuning factor is not based on G1)
- `ind` Sunlit leaf index within the leaf angular distribution

"""
function empirical_equation end;

empirical_equation(sm::BallBerrySM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;

    return G0 .+ β .* G1 .* air.s_aux.ps[3] ./ saturation_vapor_pressure(air.s_aux.t) .* leaf.flux.auxil.a_n * FT(1e-6) ./ leaf.flux.auxil.p_CO₂_s .* air.state.p_air
);
