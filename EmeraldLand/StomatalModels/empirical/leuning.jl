# This file contains the function for Leuning stomatal model

empirical_equation(sm::LeuningSM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; D0, G0, G1) = sm;

    γ = leaf.photosystem.auxil.γ_star;
    d = max(1, saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000) - air.s_aux.ps[3]);

    return G0 .+ β .* G1 ./ (1 + d / D0) .* leaf.flux.auxil.a_n .* FT(1e-6) ./ (leaf.flux.auxil.p_CO₂_s .- γ) .* air.state.p_air
);
