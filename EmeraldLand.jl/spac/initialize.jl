#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-27: add function to initialize SPAC
#     2022-Jun-27: add leaf area controller to make sure soil and leaf areas are consistent with leaf area index
#     2023-Mar-27: initialize soil and leaf e as well (because T, SWC may be changed)
#     2023-Apr-13: add config to function call
#     2023-May-19: use δlai per canopy layer
#     2023-Jun-12: initialize soil trace gas as well
#
#######################################################################################################################################################################################################
"""

    initialize!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}) where {FT<:AbstractFloat}

Initialize the SPAC, given
- `spac` `MultiLayerSPAC` SPAC
- `config` Configurations of spac model

"""
function initialize! end

initialize!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}) where {FT<:AbstractFloat} = (
    (; CANOPY, DIM_LAYER, LEAVES, SOIL) = spac;

    # make sure soil energy is correctly scaled with temperature and soil water content
    for _slayer in SOIL.LAYERS
        _slayer.e = (_slayer.CP * _slayer.ρ + _slayer.θ * CP_L() * ρ_H₂O()) * _slayer.t;
        _slayer.TRACES.n_H₂O = saturation_vapor_pressure(_slayer.t) * _slayer.ΔZ * max(0, _slayer.VC.Θ_SAT - _slayer.θ) / (GAS_R() * _slayer.t);
        _slayer.TRACES.n_N₂  = spac.AIR[1].P_AIR * 0.79 * _slayer.ΔZ * max(0, _slayer.VC.Θ_SAT - _slayer.θ) / (GAS_R() * _slayer.t);
        _slayer.TRACES.n_O₂  = spac.AIR[1].P_AIR * 0.209 * _slayer.ΔZ * max(0, _slayer.VC.Θ_SAT - _slayer.θ) / (GAS_R() * _slayer.t);
        _slayer.TRACES.p_H₂O = saturation_vapor_pressure(_slayer.t);
        _slayer.TRACES.p_N₂  = spac.AIR[1].P_AIR * 0.79;
        _slayer.TRACES.p_O₂  = spac.AIR[1].P_AIR * 0.209;
    end;

    # make sure leaf area index setup and energy are correct
    for _i in 1:DIM_LAYER
        _clayer = LEAVES[_i];
        _clayer.HS.AREA = SOIL.AREA * CANOPY.δlai[_i];
        _clayer.e = (_clayer.CP * _clayer.BIO.lma * 10 + _clayer.HS.v_storage * CP_L_MOL(FT)) * _clayer.t;
    end;

    # initialize leaf level spectra
    leaf_spectra!(spac, config);

    # initialize stomatal conductance
    stomatal_conductance!(spac, FT(0));

    return nothing
);
