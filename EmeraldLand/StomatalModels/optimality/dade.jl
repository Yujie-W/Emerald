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
#
#######################################################################################################################################################################################################
"""

    ∂A∂E(leaf::Union{CanopyLayer{FT}, Leaf{FT}}, air::AirLayer{FT}) where {FT}
    ∂A∂E(leaf::Union{CanopyLayer{FT}, Leaf{FT}}, air::AirLayer{FT}, ind::Int) where {FT}

Return the partial derivative of A per E, given
- `leaf` `Leaf` type leaf
- `air` `AirLayer` type environmental conditions
- `ind` Index of the leaf sunlit part

"""
function ∂A∂E end;

∂A∂E(leaf::CanopyLayer{FT}, air::AirLayer{FT}) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the A and E at the current setting
    gs1 = [leaf.flux.state.g_H₂O_s_sunlit[:]; leaf.flux.state.g_H₂O_s_shaded];
    gh1 = 1 ./ (1 ./ gs1 .+ 1 ./ (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e1  = gh1 * d / air.state.p_air;
    a1  = [leaf.flux.auxil.a_n_sunlit[:]; leaf.flux.auxil.a_n_shaded];

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    gs2 = gs1 .+ FT(0.0001);
    gh2 = 1 ./ (1 ./ gs2 .+ 1 ./ (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    gc2 = 1 ./ (FT(1.6) ./ gs2 .+ 1 ./ leaf.flux.auxil.g_CO₂_b);
    e2  = gh2 * d / air.state.p_air;
    a2  = photosynthesis_only!(leaf.photosystem, air, gc2, [leaf.flux.auxil.ppar_sunlit[:]; leaf.flux.auxil.ppar_shaded]);

    return (a2 - a1) / (e2 - e1)
);

∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the A and E at the current setting
    gs1 = leaf.flux.state.g_H₂O_s_shaded;
    gh1 = 1 / (1 / gs1 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e1  = gh1 * d / air.state.p_air;
    a1  = leaf.flux.auxil.a_n_shaded;

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    gs2 = gs1 + FT(0.0001);
    gh2 = 1 / (1 / gs2 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    gc2 = 1 / (FT(1.6) / gs2 + 1 / leaf.flux.auxil.g_CO₂_b);
    e2  = gh2 * d / air.state.p_air;
    a2  = photosynthesis_only!(leaf.photosystem, air, gc2, leaf.flux.auxil.ppar_shaded);

    return (a2 - a1) / (e2 - e1)
);

∂A∂E(leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the A and E at the current setting
    gs1 = leaf.flux.state.g_H₂O_s_sunlit[ind];
    gh1 = 1 / (1 / gs1 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e1  = gh1 * d / air.state.p_air;
    a1  = leaf.flux.auxil.a_n_sunlit[ind];

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    gs2 = gs1 + FT(0.0001);
    gh2 = 1 / (1 / gs2 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    gc2 = 1 / (FT(1.6) / gs2 + 1 / leaf.flux.auxil.g_CO₂_b);
    e2 = gh2 * d / air.state.p_air;
    a2 = photosynthesis_only!(leaf.photosystem, air, gc2, leaf.flux.auxil.ppar_sunlit[ind]);

    return (a2 - a1) / (e2 - e1)
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
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

    # compute the A and E at the current setting
    gs1 = [leaf.flux.state.g_H₂O_s_sunlit[:]; leaf.flux.state.g_H₂O_s_shaded];
    gh1 = 1 ./ (1 ./ gs1 .+ 1 ./ (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e1  = gh1 * d / air.state.p_air;
    a1  = [leaf.flux.auxil.a_n_sunlit[:]; leaf.flux.auxil.a_n_shaded];

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    gs2 = gs1 .+ FT(0.0001);
    gh2 = 1 ./ (1 ./ gs2 .+ 1 ./ (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    gc2 = 1 ./ (FT(1.6) ./ gs2 .+ 1 ./ leaf.flux.auxil.g_CO₂_b);
    e2  = gh2 * d / air.state.p_air;
    a2  = photosynthesis_only!(leaf.photosystem, air, gc2, [leaf.flux.auxil.ppar_sunlit[:]; leaf.flux.auxil.ppar_shaded]);

    dade = (a2 .- a1) ./ (e2 .- e1);

    leaf.flux.auxil.∂A∂E_sunlit[:] .= dade[1:end-1];
    leaf.flux.auxil.∂A∂E_shaded = dade[end];

    return nothing
);

∂A∂E!(cache::SPACCache{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
    p_s = saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000);
    d = max(1, p_s - air.s_aux.ps[3]);

    # compute the A and E at the current setting for shaded leaves
    gs1 = leaf.flux.state.g_H₂O_s_shaded;
    gh1 = 1 / (1 / gs1 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e1  = gh1 * d / air.state.p_air;
    a1  = leaf.flux.auxil.a_n_shaded;

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    gs2 = gs1 + FT(0.0001);
    gh2 = 1 / (1 / gs2 + 1 / (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    gc2 = 1 / (FT(1.6) / gs2 + 1 / leaf.flux.auxil.g_CO₂_b);
    e2  = gh2 * d / air.state.p_air;
    a2  = photosynthesis_only!(leaf.photosystem, air, gc2, leaf.flux.auxil.ppar_shaded);

    leaf.flux.auxil.∂A∂E_shaded = (a2 - a1) / (e2 - e1);

    # compute the A and E at the current setting for sunlit leaves
    gs1  = leaf.flux.state.g_H₂O_s_sunlit;
    gh1  = cache.cache_incl_azi_1;
    e1   = cache.cache_incl_azi_2;
    gh1 .= 1 ./ (1 ./ gs1 .+ 1 ./ (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    e1  .= gh1 .* (d / air.state.p_air);
    a1   = leaf.flux.auxil.a_n_sunlit;

    # compute the A and E when g_sw increases by 0.0001 mol m⁻² s⁻¹
    gs2  = cache.cache_incl_azi_3;
    gh2  = cache.cache_incl_azi_4;
    gc2  = cache.cache_incl_azi_5;
    e2   = cache.cache_incl_azi_6;
    a2   = cache.cache_incl_azi_7;
    gs2 .= gs1 .+ FT(0.0001);
    gh2 .= 1 ./ (1 ./ gs2 .+ 1 ./ (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    gc2 .= 1 ./ (FT(1.6) ./ gs2 .+ 1 ./ leaf.flux.auxil.g_CO₂_b);
    e2  .= gh2 .* (d / air.state.p_air);
    a2  .= photosynthesis_only!.((leaf.photosystem,), (air,), gc2, leaf.flux.auxil.ppar_sunlit); # need to run the shaded leaves first to update temperature dependent variables

    leaf.flux.auxil.∂A∂E_sunlit .= (a2 .- a1) ./ (e2 .- e1);

    return nothing
);
