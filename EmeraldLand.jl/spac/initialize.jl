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
#     2023-Jun-13: update N₂ and O₂ based on soil water content
#     2023-Jun-13: add soil gas energy into soil e
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Sep-07: add ALLOW_SOIL_EVAPORATION check
#     2023-Sep-11: rename ALLOW_SOIL_EVAPORATION to ENABLE_SOIL_EVAPORATION
#
#######################################################################################################################################################################################################
"""

    initialize!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Initialize the SPAC, given
- `config` Configurations of spac model
- `spac` `MultiLayerSPAC` SPAC

"""
function initialize! end

initialize!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; ENABLE_SOIL_EVAPORATION) = config;
    (; CANOPY, LEAVES, SOIL) = spac;

    # make sure soil energy is correctly scaled with temperature and soil water content
    for _slayer in SOIL.LAYERS
        if ENABLE_SOIL_EVAPORATION
            _δθ = max(0, _slayer.VC.Θ_SAT - _slayer.θ);
            _rt = GAS_R(FT) * _slayer.t;
            _slayer.TRACES.n_H₂O = saturation_vapor_pressure(_slayer.t, _slayer.ψ * 1000000) * _slayer.ΔZ * _δθ / _rt;
            _slayer.TRACES.n_N₂  = spac.AIR[1].P_AIR * 0.79 * _slayer.ΔZ * _δθ / _rt;
            _slayer.TRACES.n_O₂  = spac.AIR[1].P_AIR * 0.209 * _slayer.ΔZ * _δθ / _rt;
        end;
        _cp_gas = (_slayer.TRACES.n_H₂O * CP_V_MOL(FT) + (_slayer.TRACES.n_CH₄ + _slayer.TRACES.n_CO₂ + _slayer.TRACES.n_N₂ + _slayer.TRACES.n_O₂) * CP_D_MOL(FT)) / _slayer.ΔZ;
        _slayer.e = (_slayer.ρ * _slayer.CP + _slayer.θ * ρ_H₂O(FT) * CP_L(FT) + _cp_gas) * _slayer.t;
    end;

    # make sure leaf area index setup and energy are correct
    for _i in eachindex(LEAVES)
        _clayer = LEAVES[_i];
        _clayer.HS.AREA = SOIL.AREA * CANOPY.δlai[_i];
        _clayer.e = (_clayer.CP * _clayer.BIO.state.lma * 10 + _clayer.HS.v_storage * CP_L_MOL(FT)) * _clayer.t;
    end;

    # initialize leaf level spectra
    leaf_spectra!(config, spac);

    # initialize stomatal conductance
    stomatal_conductance!(spac, FT(0));

    return nothing
);
