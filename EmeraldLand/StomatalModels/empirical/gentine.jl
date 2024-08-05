# This file contains the function for Gentine stomatal model

empirical_equation(sm::GentineSM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;

    return G0 .+ β .* G1 .* leaf.flux.auxil.a_n .* FT(1e-6) ./ leaf.flux.auxil.p_CO₂_i .* air.state.p_air
);
