# This file contains functions to compute the partial derivatives of A versus E

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: migrate function from older version
#     2022-Jul-07: add the method for Leaf
#     2022-Jul-07: add the method for Leaf (shaded leaf)
#     2022-Jul-07: add the method for Leaf (sunlit leaf)
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Aug-23: add nan check for some methods (not all, to do later)
#     2023-Aug-27: make sure D >= 1 when compute dAdE
# To do
#     TODO: figure out why the ratios are 1.35 and 1.6, and make them more accurate
#
#######################################################################################################################################################################################################
"""

    ∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT}
    ∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int) where {FT}

Return the partial derivative of A per E, given
- `leaf` `Leaf` type leaf
- `air` `AirLayer` type environmental conditions
- `ind` Index of the leaf sunlit part

"""
function ∂A∂E end;

∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    d = max(1, p_s - air.auxil.ps[3]);

    # compute the A and E at the current setting
    gs1 = leaf.flux.state.g_H₂O_s_shaded;
    gh1 = 1 / (1 / gs1 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e1  = gh1 * d / air.state.p_air;
    a1  = leaf.flux.auxil.a_n_shaded;

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    gs2 = gs1 + FT(0.0001);
    gh2 = 1 / (1 / gs2 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    gc2 = 1 / (FT(1.6) / gs2 + 1 / leaf.flux.auxil.g_CO₂_b);
    photosynthesis_only!(leaf.photosystem, air, gc2, leaf.flux.auxil.ppar_shaded, leaf.energy.auxil.t);
    e2 = gh2 * d / air.state.p_air;
    a2 = leaf.photosystem.auxil.a_n;

    return (a2 - a1) / (e2 - e1)
);

∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000);
    d = max(1, p_s - air.auxil.ps[3]);

    # compute the A and E at the current setting
    gs1 = leaf.flux.state.g_H₂O_s_sunlit[ind];
    gh1 = 1 / (1 / gs1 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e1  = gh1 * d / air.state.p_air;
    a1  = leaf.flux.auxil.a_n_sunlit[ind];

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    gs2 = gs1 + FT(0.0001);
    gh2 = 1 / (1 / gs2 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    gc2 = 1 / (FT(1.6) / gs2 + 1 / leaf.flux.auxil.g_CO₂_b);
    photosynthesis_only!(leaf.photosystem, air, gc2, leaf.flux.auxil.ppar_sunlit[ind], leaf.energy.auxil.t);
    e2 = gh2 * d / air.state.p_air;
    a2 = leaf.photosystem.auxil.a_n;

    return (a2 - a1) / (e2 - e1);
);
