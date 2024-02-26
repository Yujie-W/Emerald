# This file contains the function for Medlyn stomatal model

empirical_equation(sm::MedlynSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;

    d = max(1, saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000) - air.s_aux.ps[3]);

    return G0 + FT(1.6) * (1 + β * G1 / sqrt(d)) * leaf.flux.auxil.a_n_shaded * FT(1e-6) / air.s_aux.ps[2] * air.state.p_air
);

empirical_equation(sm::MedlynSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}, ind::Int; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;

    d = max(1, saturation_vapor_pressure(leaf.energy.auxil.t, leaf.capacitor.auxil.p_leaf * 1000000) - air.s_aux.ps[3]);

    return G0 + FT(1.6) * (1 + β * G1 / sqrt(d)) * leaf.flux.auxil.a_n_sunlit[ind] * FT(1e-6) / air.s_aux.ps[2] * air.state.p_air
);
