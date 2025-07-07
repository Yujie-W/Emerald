# This file contains the function for Medlyn stomatal model

empirical_equation(sm::MedlynSM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    (; G0, G1) = sm;

    d = max(1, saturation_vapor_pressure(leaf.energy.s_aux.t, leaf.capacitor.state.p_leaf * 1000000) - air.s_aux.ps[3]);

    return (G0 .+ FT(1.6) .* (1 + G1 / sqrt(d)) .* leaf.flux.auxil.a_n * FT(1e-6) ./ air.s_aux.ps[2] .* air.state.p_air) .* β
);
