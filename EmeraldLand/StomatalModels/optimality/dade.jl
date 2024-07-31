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
#     2024-Jul-24: add function to update the ∂A∂E for sunlit leaves (matrix version)
#
#######################################################################################################################################################################################################
"""

    ∂A∂E!(cache::SPACCache{FT}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, air::AirLayer{FT}) where {FT}

Update the ∂A∂E for sunlit leaves, given
- `cache` `SPACCache` type cache
- `leaf` `Leaf` type leaf
- `air` `AirLayer` type environmental conditions

"""
function ∂A∂E! end;

∂A∂E!(cache::SPACCache{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # unpack the cache vars (note do not reuse these vars in the photosynthesis_only! function)
    gh1 = cache.cache_incl_azi_1_1;
    e1  = cache.cache_incl_azi_1_2;
    gs2 = cache.cache_incl_azi_1_3;
    gh2 = cache.cache_incl_azi_1_4;
    gc2 = cache.cache_incl_azi_1_5;
    e2  = cache.cache_incl_azi_1_6;

    # compute the A and E at the current setting
    gs1 = leaf.flux.state.g_H₂O_s;
    gh1 .= 1 ./ (1 ./ gs1 .+ 1 ./ (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e1  .= gh1 .* d ./ air.state.p_air;
    a1  = leaf.flux.auxil.a_n;

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    gs2 .= gs1 .+ FT(0.0001);
    gh2 .= 1 ./ (1 ./ gs2 .+ 1 ./ (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    gc2 .= 1 ./ (FT(1.6) ./ gs2 .+ 1 ./ leaf.flux.auxil.g_CO₂_b);
    e2  .= gh2 .* d ./ air.state.p_air;
    a2   = photosynthesis_only!(cache, leaf.photosystem, air, gc2, leaf.flux.auxil.ppar);

    leaf.flux.auxil.∂A∂E .= (a2 .- a1) ./ (e2 .- e1);

    return nothing
);
