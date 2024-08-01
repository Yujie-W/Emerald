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
#     2024-Jul-27: set η_c and η_l to min(η_c, η_l) to avoid negative values of η
#     2024-Aug-01: generalize the function for GeneralC3Trait and GeneralC4Trait
#
#######################################################################################################################################################################################################
"""

    photosystem_temperature_dependence!(config::SPACConfiguration{FT}, ps::LeafPhotosystem{FT}, air::AirLayer{FT}, t::FT) where {FT}

Update the temperature dependencies of C3 photosynthesis model, given
- `config` `SPACConfiguration` structure
- `psm` `LeafPhotosystem` or `CanopyLayerPhotosystem` structure
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `t` Target temperature in `K`

"""
function photosystem_temperature_dependence! end;

photosystem_temperature_dependence!(
            config::SPACConfiguration{FT},
            ps::Union{CanopyLayerPhotosystem{FT}, LeafPhotosystem{FT}},
            air::AirLayer{FT},
            t::FT) where {FT} = photosystem_temperature_dependence!(config, ps.trait, ps.auxil, air, t);

photosystem_temperature_dependence!(
            config::SPACConfiguration{FT},
            pst::GeneralC3Trait{FT},
            psa::Union{CanopyLayerPhotosystemAuxil{FT}, LeafPhotosystemAuxil{FT}},
            air::AirLayer{FT},
            t::FT) where {FT} = photosystem_temperature_dependence!(config, pst, psa, pst.ACM, pst.AJM, pst.APM, air, t);

photosystem_temperature_dependence!(
            config::SPACConfiguration{FT},
            pst::GeneralC3Trait{FT},
            psa::Union{CanopyLayerPhotosystemAuxil{FT}, LeafPhotosystemAuxil{FT}},
            acm::AcMethodC3VcmaxPi,
            ajm::AjMethodC3JmaxPi,
            apm::Union{ApMethodC3Inf, ApMethodC3Vcmax},
            air::AirLayer{FT},
            t::FT) where {FT} = (
    (; ENABLE_KD_TD, PS_RATE_CONSTANTS) = config;

    psa.j_max  = pst.j_max25 * temperature_correction(pst.TD_JMAX, t);
    psa.k_c    = temperature_corrected_value(pst.TD_KC, t);
    psa.k_o    = temperature_corrected_value(pst.TD_KO, t);
    psa.r_d    = pst.r_d25 * temperature_correction(pst.TD_R, t);
    psa.v_cmax = pst.v_cmax25 * temperature_correction(pst.TD_VCMAX, t);
    psa.γ_star = temperature_corrected_value(pst.TD_Γ, t);
    psa.k_m    = psa.k_c * (1 + air.state.p_air * F_O₂(FT) / psa.k_o);

    # TODO: make this a method?
    psa.k_d = ENABLE_KD_TD ? max(0.8738, 0.0301 * (t - 273.15) + 0.0773) : PS_RATE_CONSTANTS.K_D;
    psa.ϕ_psii_max = PS_RATE_CONSTANTS.K_P / (psa.k_d + PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P);

    return nothing
);

photosystem_temperature_dependence!(
            config::SPACConfiguration{FT},
            pst::GeneralC3Trait{FT},
            psa::Union{CanopyLayerPhotosystemAuxil{FT}, LeafPhotosystemAuxil{FT}},
            acm::AcMethodC3VcmaxPi,
            ajm::AjMethodC3VqmaxPi,
            apm::Union{ApMethodC3Inf, ApMethodC3Vcmax},
            air::AirLayer{FT},
            t::FT) where {FT} = (
    (; FIX_ETA_TD, PSI_RATE_CONSTANTS) = config;

    psa.k_c    = temperature_corrected_value(pst.TD_KC, t);
    psa.k_o    = temperature_corrected_value(pst.TD_KO, t);
    psa.k_q    = temperature_corrected_value(pst.TD_KQ, t);
    psa.r_d    = pst.r_d25 * temperature_correction(pst.TD_R, t);
    psa.v_cmax = pst.v_cmax25 * temperature_correction(pst.TD_VCMAX, t);
    psa.γ_star = temperature_corrected_value(pst.TD_Γ, t);
    psa.η_c    = temperature_corrected_value(pst.TD_ηC, t);
    psa.η_l    = temperature_corrected_value(pst.TD_ηL, t);
    psa.k_m    = psa.k_c * (1 + air.state.p_air * F_O₂(FT) / psa.k_o);
    psa.v_qmax = pst.b₆f * psa.k_q;

    if FIX_ETA_TD
        psa.η_c = min(pst.TD_ηC.VAL_REF, psa.η_c);
        psa.η_l = min(pst.TD_ηL.VAL_REF, psa.η_c);
    end;

    psa.ϕ_psi_max = PSI_RATE_CONSTANTS.K_P / (PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_P);

    return nothing
);

photosystem_temperature_dependence!(
            config::SPACConfiguration{FT},
            pst::GeneralC3Trait{FT},
            psa::Union{CanopyLayerPhotosystemAuxil{FT}, LeafPhotosystemAuxil{FT}},
            acm::AcMethodC4Vcmax,
            ajm::AjMethodC4JPSII,
            apm::ApMethodC4VcmaxPi,
            air::AirLayer{FT},
            t::FT) where {FT} = (
    (; ENABLE_KD_TD, PS_RATE_CONSTANTS) = config;

    psa.k_pep_clm = temperature_corrected_value(pst.TD_KPEP, t);
    psa.r_d       = pst.r_d25 * temperature_correction(pst.TD_R, t);
    psa.v_cmax    = pst.v_cmax25 * temperature_correction(pst.TD_VCMAX, t);

    # TODO: make this a method?
    psa.k_d = ENABLE_KD_TD ? max(0.8738, 0.0301 * (t - 273.15) + 0.0773) : PS_RATE_CONSTANTS.K_D;
    psa.ϕ_psii_max = PS_RATE_CONSTANTS.K_P / (psa.k_d + PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P);

    return nothing
);

photosystem_temperature_dependence!(
            config::SPACConfiguration{FT},
            pst::GeneralC3Trait{FT},
            psa::Union{CanopyLayerPhotosystemAuxil{FT}, LeafPhotosystemAuxil{FT}},
            acm::AcMethodC4Vcmax,
            ajm::AjMethodC4JPSII,
            apm::ApMethodC4VpmaxPi,
            air::AirLayer{FT},
            t::FT) where {FT} = (
    (; ENABLE_KD_TD, PS_RATE_CONSTANTS) = config;

    psa.k_pep  = temperature_corrected_value(pst.TD_KPEP, t);
    psa.r_d    = pst.r_d25 * temperature_correction(pst.TD_R, t);
    psa.v_cmax = pst.v_cmax25 * temperature_correction(pst.TD_VCMAX, t);
    psa.v_pmax = pst.v_pmax25 * temperature_correction(pst.TD_VPMAX, t);

    # TODO: make this a method?
    psa.k_d = ENABLE_KD_TD ? max(0.8738, 0.0301 * (t - 273.15) + 0.0773) : PS_RATE_CONSTANTS.K_D;
    psa.ϕ_psii_max = PS_RATE_CONSTANTS.K_P / (psa.k_d + PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P);

    return nothing
);
