#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Jul-19: add functions to obtain RMSE of A-Ci curve
#     2024-Jul-22: add support to C3CLM, C3FvCB, and C3VJP
#     2024-Jul-27: add option to turn on/off Rd fitting
#     2024-Jul-27: fit Γ_star as well for C3 models
#     2024-Aug-01: use GeneralC3Trait and GeneralC4Trait
#
#######################################################################################################################################################################################################
"""

    aci_rmse(config::SPACConfiguration{FT},
             ps::LeafPhotosystem{FT},
             air::AirLayer{FT},
             df::DataFrame,
             params::Vector;
             fit_rd::Bool = false,
             fitted_γ::Union{Nothing, Number} = nothing) where {FT}

Compute the RMSE of A-Ci curve (will be abstractized using the trait and methods embedded), given
- `config` `SPACConfiguration` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `params` Vector of parameters (Vcmax25, Vpmax25, b₆f, Jmax25, Rd25 depending on the photosynthesis model)
- `fit_rd` Fit Rd or not (last fitting parameter; default: false)
- `fitted_γ` Fitted Γ* value, will be prescribed if given (default: nothing)

"""
function aci_rmse end;

aci_rmse(config::SPACConfiguration{FT},
         ps::LeafPhotosystem{FT},
         air::AirLayer{FT},
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Union{Nothing, Number} = nothing) where {FT} = aci_rmse(config, ps, ps.trait, air, df, params; fit_rd = fit_rd, fitted_γ = fitted_γ);

aci_rmse(config::SPACConfiguration{FT},
         ps::LeafPhotosystem{FT},
         pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
         air::AirLayer{FT},
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Union{Nothing, Number} = nothing) where {FT} = aci_rmse(config, ps, pst, pst.ACM, pst.AJM, pst.APM, air, df, params; fit_rd = fit_rd, fitted_γ = fitted_γ);

aci_rmse(config::SPACConfiguration{FT},
         ps::LeafPhotosystem{FT},
         pst::GeneralC3Trait{FT},
         acm::AcMethodC3VcmaxPi,
         ajm::AjMethodC3JmaxPi,
         apm::ApMethodC3Vcmax,
         air::AirLayer{FT},
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Union{Nothing, Number} = nothing) where {FT} = (
    pst.v_cmax25 = params[1];
    pst.j_max25 = params[2];
    if isnothing(fitted_γ)
        pst.TD_Γ.VAL_REF = params[3];
    else
        pst.TD_Γ.VAL_REF = fitted_γ;
    end;
    pst.r_d25 = fit_rd ? params[4] : params[1] * 0.015;

    return rmse(aci_curve(config, ps, air, df), df.A_NET)
);

aci_rmse(config::SPACConfiguration{FT},
         ps::LeafPhotosystem{FT},
         pst::GeneralC3Trait{FT},
         acm::AcMethodC3VcmaxPi,
         ajm::AjMethodC3VqmaxPi,
         apm::ApMethodC3Vcmax,
         air::AirLayer{FT},
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Union{Nothing, Number} = nothing) where {FT} = (
    pst.v_cmax25 = params[1];
    pst.b₆f = params[2];
    if isnothing(fitted_γ)
        pst.TD_Γ.VAL_REF = params[3];
    else
        pst.TD_Γ.VAL_REF = fitted_γ;
    end;
    pst.r_d25 = fit_rd ? params[4] : params[1] * 0.015;

    return rmse(aci_curve(config, ps, air, df), df.A_NET)
);

aci_rmse(config::SPACConfiguration{FT},
         ps::LeafPhotosystem{FT},
         pst::GeneralC4Trait{FT},
         acm::AcMethodC4Vcmax,
         ajm::AjMethodC4JPSII,
         apm::ApMethodC4VcmaxPi,
         air::AirLayer{FT},
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Nothing = nothing) where {FT} = (
    pst.v_cmax25 = params[1];
    pst.r_d25 = fit_rd ? params[2] : params[1] * 0.015;

    return rmse(aci_curve(config, ps, air, df), df.A_NET)
);

aci_rmse(config::SPACConfiguration{FT},
         ps::LeafPhotosystem{FT},
         pst::GeneralC4Trait{FT},
         acm::AcMethodC4Vcmax,
         ajm::AjMethodC4JPSII,
         apm::ApMethodC4VpmaxPi,
         air::AirLayer{FT},
         df::DataFrame,
         params::Vector;
         fit_rd::Bool = false,
         fitted_γ::Nothing = nothing) where {FT} = (
    pst.v_cmax25 = params[1];
    pst.v_pmax25 = params[2];
    pst.r_d25 = fit_rd ? params[3] : params[1] * 0.015;

    return rmse(aci_curve(config, ps, air, df), df.A_NET)
);
