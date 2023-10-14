# This file contains the function for Leuning stomatal model

empirical_equation(sm::LeuningSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; D0, G0, G1) = sm;

    γ = (typeof(leaf.photosystem) <: C4VJP) ? 0 : leaf.photosystem.auxil.γ_star;
    d = max(1, saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000) - air.auxil.ps[3]);

    return G0 + β * G1 / (1 + d / D0) * leaf.flux.auxil.a_n_shaded * FT(1e-6) / (leaf.flux.auxil.p_CO₂_s_shaded- γ) * air.state.p_air
);

empirical_equation(sm::LeuningSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT} = (
    (; D0, G0, G1) = sm;

    γ = (typeof(leaf.photosystem) <: C4VJP) ? 0 : leaf.photosystem.auxil.γ_star;
    d = max(1, saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000) - air.auxil.ps[3]);

    return G0 + β * G1 / (1 + d / D0) * leaf.flux.auxil.a_n_sunlit[ind] * FT(1e-6) / (leaf.flux.auxil.p_CO₂_s_sunlit[ind] - γ) * air.state.p_air
);
