using Emerald
using Emerald.Namespace
using Emerald.SPAC


FT = Float64;
config = Emerald.Namespace.SPACConfig(FT);
        config.METHODS.C3_AC_METHOD = Namespace.AcMethodC3VcmaxPi();
        config.METHODS.C3_AJ_METHOD = Namespace.AjMethodC3VqmaxPi();
        config.METHODS.C3_AP_METHOD = Namespace.ApMethodC3Vcmax();
        config.METHODS.COLIMIT_J = Namespace.SerialColimit();
        config.METHODS.TD_VCMAX_C3.ΔHA = 63000;
        config.METHODS.TD_VCMAX_C3.ΔHD = 204000;
        config.METHODS.TD_KQ = Namespace.ArrheniusPeak{Float64}(298.15, 300, 28500, 223500, 700);
        config.METHODS.TD_Γ.VAL_REF = 4.56;
        config.METHODS.TD_Γ.ΔHA = 11800;
        config.METHODS.TD_ηC.VAL_REF = 4 * 3 / 14;
        config.METHODS.TD_ηC.ΔHA = 28500;
        config.METHODS.TD_ηC.ΔHD = 223500;
        config.METHODS.TD_ηC.ΔSV = 700;
        config.METHODS.TD_ηL.VAL_REF = 3 * 3 / 14;
        config.METHODS.TD_ηL.ΔHA = 28500;
        config.METHODS.TD_ηL.ΔHD = 223500;
        config.METHODS.TD_ηL.ΔSV = 700;
config.METHODS.FLUORESCENCE_METHOD_C3 = Namespace.CytochromeFluorescenceModel();

spac = Namespace.BulkSPAC(config);
SPAC.initialize_spac!(config, spac);
SPAC.spac!(config, spac, FT(1));
