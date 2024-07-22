# This file contains function to make temperature dependencies for photosynthesis

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-13: use ClimaCache types, which uses ΔHA, ΔHD, and ΔSV directly
#     2022-Jan-13: add optional input t_ref to allow for manually setting reference temperature
#     2022-Jul-29: add support to Q10Peak
#
#######################################################################################################################################################################################################
"""

    temperature_correction(td::Arrhenius{FT}, t::FT; t_ref::FT = td.T_REF) where {FT}
    temperature_correction(td::ArrheniusPeak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT}
    temperature_correction(td::Q10{FT}, t::FT; t_ref::FT = td.T_REF) where {FT}
    temperature_correction(td::Q10Peak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT}

Return the correction ratio for a temperature dependent variable, given
- `td` `Arrhenius`, `ArrheniusPeak`, `Q10`, or `Q10Peak` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)

"""
function temperature_correction end;

temperature_correction(td::Arrhenius{FT}, t::FT; t_ref::FT = td.T_REF) where {FT} = exp( td.ΔHA / GAS_R(FT) * (1/t_ref - 1/t) );

temperature_correction(td::ArrheniusPeak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT} = (
    (; ΔHA, ΔHD, ΔSV) = td;

    # f_a: activation correction, f_b: de-activation correction
    f_a = exp( ΔHA / GAS_R(FT) * (1 / t_ref - 1 / t) );
    f_b = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t_ref))) / (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t)));

    return f_a * f_b
);

temperature_correction(td::Q10{FT}, t::FT; t_ref::FT = td.T_REF) where {FT} = td.Q_10 ^ ( (t - t_ref) / 10 );

temperature_correction(td::Q10Peak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT} = (
    (; ΔHD, ΔSV) = td;

    # f_a: activation correction, f_b: de-activation correction
    f_a = td.Q_10 ^ ( (t - t_ref) / 10 );
    f_b = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t_ref))) / (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t)));

    return f_a * f_b
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-13: use ClimaCache types, which uses ΔHA, ΔHD, and ΔSV directly
#     2022-Jan-13: add optional input t_ref to allow for manually setting reference temperature
#     2022-Jul-29: add support to Q10Peak
#
#######################################################################################################################################################################################################
"""

    temperature_corrected_value(td::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}}, t::FT; t_ref::FT = td.T_REF) where {FT}

Return the temperature corrected value, given
- `td` `Arrhenius`, `ArrheniusPeak`, `Q10`, or `Q10Peak` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)

"""
function temperature_corrected_value(td::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}}, t::FT; t_ref::FT = td.T_REF) where {FT}
    return td.VAL_REF * temperature_correction(td, t; t_ref=t_ref)
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: use ClimaCache types, which saves photosystem temperature dependencies within the struct
#     2022-Feb-07: add method for C3CytochromeModel photosynthesis model
#     2022-Feb-07: add v_qmax without temperature dependency
#     2022-Mar-01: add temperature dependencies for k_q, v_qmax, η_c, and η_l
#     2024-Apr-15: add support for C4CLMTrait model
#
#######################################################################################################################################################################################################
"""

    photosystem_temperature_dependence!(psm::C3Cyto{FT}, air::AirLayer{FT}, t::FT) where {FT}
    photosystem_temperature_dependence!(psm::C3VJP{FT}, air::AirLayer{FT}, t::FT) where {FT}
    photosystem_temperature_dependence!(psm::C4VJP{FT}, air::AirLayer{FT}, t::FT) where {FT}

Update the temperature dependencies of C3 photosynthesis model, given
- `psm` `C3Cyto`, `C3VJP`, or `C4VJP` type photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `t` Target temperature in `K`

"""
function photosystem_temperature_dependence! end;

photosystem_temperature_dependence!(psm::LeafPhotosystem{FT}, air::AirLayer{FT}, t::FT) where {FT} = photosystem_temperature_dependence!(psm.trait, psm.auxil, air, t);

photosystem_temperature_dependence!(pst::Union{C3CytoTrait{FT}, C3JBTrait{FT}}, psa::LeafPhotosystemAuxil{FT}, air::AirLayer{FT}, t::FT) where {FT} = (
    psa.k_c    = temperature_corrected_value(pst.TD_KC, t);
    psa.k_o    = temperature_corrected_value(pst.TD_KO, t);
    psa.k_q    = temperature_corrected_value(pst.TD_KQ, t);
    psa.γ_star = temperature_corrected_value(pst.TD_Γ , t);
    psa.η_c    = temperature_corrected_value(pst.TD_ηC, t);
    psa.η_l    = temperature_corrected_value(pst.TD_ηL, t);
    psa.r_d    = pst.r_d25    * temperature_correction(pst.TD_R, t);
    psa.v_cmax = pst.v_cmax25 * temperature_correction(pst.TD_VCMAX, t);
    psa.k_m    = psa.k_c * (1 + air.state.p_air * F_O₂(FT) / psa.k_o);
    psa.v_qmax = pst.b₆f * psa.k_q;

    psa.ϕ_psi_max = pst.K_PSI / (pst.K_D + pst.K_F + pst.K_PSI);

    return nothing
);

photosystem_temperature_dependence!(pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}}, psa::LeafPhotosystemAuxil{FT}, air::AirLayer{FT}, t::FT) where {FT} = (
    psa.k_c    = temperature_corrected_value(pst.TD_KC, t);
    psa.k_o    = temperature_corrected_value(pst.TD_KO, t);
    psa.γ_star = temperature_corrected_value(pst.TD_Γ , t);
    psa.r_d    = pst.r_d25    * temperature_correction(pst.TD_R, t);
    psa.v_cmax = pst.v_cmax25 * temperature_correction(pst.TD_VCMAX, t);
    psa.j_max  = pst.j_max25  * temperature_correction(pst.TD_JMAX, t);
    psa.k_m    = psa.k_c * (1 + air.state.p_air * F_O₂(FT) / psa.k_o);

    # TODO: add a TD_KD in the model in the future like psd.TD_KC
    psa.k_d = max(0.8738, 0.0301 * (t - 273.15) + 0.0773);
    psa.ϕ_psii_max = pst.K_PSII / (psa.k_d + pst.K_F + pst.K_PSII);

    return nothing
);

photosystem_temperature_dependence!(pst::C4CLMTrait{FT}, psa::LeafPhotosystemAuxil{FT}, air::AirLayer{FT}, t::FT) where {FT} = (
    psa.k_pep_clm = temperature_corrected_value(pst.TD_KPEP, t);
    psa.r_d       = pst.r_d25    * temperature_correction(pst.TD_R, t);
    psa.v_cmax    = pst.v_cmax25 * temperature_correction(pst.TD_VCMAX, t);

    # TODO: add a TD_KD in the model in the future like psd.TD_KC
    psa.k_d = max(0.8738, 0.0301 * (t - 273.15) + 0.0773);
    psa.ϕ_psii_max = pst.K_PSII / (psa.k_d + pst.K_F + pst.K_PSII);

    return nothing
);

photosystem_temperature_dependence!(pst::C4VJPTrait{FT}, psa::LeafPhotosystemAuxil{FT}, air::AirLayer{FT}, t::FT) where {FT} = (
    psa.k_pep  = temperature_corrected_value(pst.TD_KPEP, t);
    psa.r_d    = pst.r_d25    * temperature_correction(pst.TD_R, t);
    psa.v_cmax = pst.v_cmax25 * temperature_correction(pst.TD_VCMAX, t);
    psa.v_pmax = pst.v_pmax25 * temperature_correction(pst.TD_VPMAX, t);

    # TODO: add a TD_KD in the model in the future like psd.TD_KC
    psa.k_d = max(0.8738, 0.0301 * (t - 273.15) + 0.0773);
    psa.ϕ_psii_max = pst.K_PSII / (psa.k_d + pst.K_F + pst.K_PSII);

    return nothing
);
