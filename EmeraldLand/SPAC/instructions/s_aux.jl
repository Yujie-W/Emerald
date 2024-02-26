#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-26: add s_aux! method for AirLayerState-dependent variables
#     2024-Feb-26: add s_aux! method for AirLayer
#
#######################################################################################################################################################################################################
s_aux!(t_aux::AirLayerTDAuxil{FT}, state::AirLayerState{FT}, s_aux::AirLayerSDAuxil{FT}) where {FT} = (
    s_aux.t = state.Σe / heat_capacitance(state);
    for i in 1:5
        s_aux.ps[i] = (state.ns[i] * GAS_R(FT) * s_aux.t) / t_aux.δz;
    end;
    s_aux.f_CO₂ = s_aux.ps[2] / state.p_air * 1e6;

    return nothing
);

s_aux!(air::AirLayer{FT}) where {FT} = s_aux!(air.t_aux, air.state, air.s_aux);
