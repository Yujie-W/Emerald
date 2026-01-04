"""

    aci_rmse(config::SPACConfig{FT},
             ps::LeafPhotosystem{FT},
             air::AirLayer{FT},
             df::DataFrame,
             prams::Vector{String},
             xxx::Vector) where {FT}

Compute the RMSE of A-Ci curve (will be abstractized using the trait and methods embedded), given
- `config` `SPACConfig` struct
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `df` DataFrame with columns `P_I`, `PPAR`, `T_LEAF`, and `A_NET`
- `params` Vector of fitting parameters (e.g., ["Vcmax25", "Jmax25", "Γstar25", "Rd25", "b₆f"])
- `xxx` Vector of parameters (Vcmax25, Vpmax25, b₆f, Jmax25, Rd25 depending on the photosynthesis model)

"""
function aci_rmse end;

aci_rmse(config::SPACConfig{FT},
         ps::LeafPhotosystem{FT},
         air::AirLayer{FT},
         df::DataFrame,
         params::Vector{String},
         xxx::Vector) where {FT} = aci_rmse(config, ps, ps.trait, air, df, params, xxx);

aci_rmse(config::SPACConfig{FT},
         ps::LeafPhotosystem{FT},
         pst::GeneralC3Trait{FT},
         air::AirLayer{FT},
         df::DataFrame,
         params::Vector{String},
         xxx::Vector) where {FT} =
    aci_rmse(config, ps, pst, config.METHODS.C3_AC_METHOD, config.METHODS.C3_AJ_METHOD, config.METHODS.C3_AP_METHOD, air, df, params, xxx);

aci_rmse(config::SPACConfig{FT},
         ps::LeafPhotosystem{FT},
         pst::GeneralC4Trait{FT},
         air::AirLayer{FT},
         df::DataFrame,
         params::Vector{String},
         xxx::Vector) where {FT} =
    aci_rmse(config, ps, pst, config.METHODS.C4_AC_METHOD, config.METHODS.C4_AJ_METHOD, config.METHODS.C4_AP_METHOD, air, df, params, xxx);
