"""
Method configuration for the SPAC model
"""
Base.@kwdef mutable struct SPACMethods{FT<:AbstractFloat}
    # fluorescence methods
    "Fluorescence method"
    FLUORESCENCE_METHOD::AbstractFluorescenceMethod{FT} = KNFluorescenceModel{FT}()
    "Fluorescence spectra method"
    FLUORESCENCE_SPECTRA_METHOD::AbstractFluorescenceSpectraMethod = PlatespectFluorescenceSpectra()

    # photosynthesis models - Ac, Aj, Ap methods
    "C3 Ac method"
    C3_AC_METHOD::AbstractAcMethod = AcMethodC3VcmaxPi()
    "C3 Aj method"
    C3_AJ_METHOD::AbstractAjMethod = AjMethodC3JmaxPi()
    "C3 Ap method"
    C3_AP_METHOD::AbstractApMethod = ApMethodC3Vcmax()
    "C4 Ac method"
    C4_AC_METHOD::AbstractAcMethod = AcMethodC4Vcmax()
    "C4 Aj method"
    C4_AJ_METHOD::AbstractAjMethod = AjMethodC4JPSII()
    "C4 Ap method"
    C4_AP_METHOD::AbstractApMethod = ApMethodC4VcmaxPi()

    # photosynthesis models - colimitation methods
    "Colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::AbstractColimit{FT} = MinimumColimit{FT}()
    "Colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::AbstractColimit{FT} = MinimumColimit{FT}()
    "Ccolimitation method for J (for C3 only)"
    COLIMIT_J::AbstractColimit{FT} = ColimitJCLM(FT)

    # photosynthesis models - temperature dependency methods
    "Jmax temperature dependency"
    TD_JMAX::AbstractTemperatureDependency{FT} = JmaxTDCLM(FT)
    "Kc temperature dependency"
    TD_KC::AbstractTemperatureDependency{FT} = KcTDCLM(FT)
    "Ko temperature dependency"
    TD_KO::AbstractTemperatureDependency{FT} = KoTDCLM(FT)
    "Kpep temperature dependency"
    TD_KPEP::AbstractTemperatureDependency{FT} = KpepTDBoyd(FT)
    "Kpep temperature dependency to use with C4CLM method"
    TD_KPEP_CLM::AbstractTemperatureDependency{FT} = Q10TDKpepCLM(FT)
    "Kq temperature dependency"
    TD_KQ::AbstractTemperatureDependency{FT} = KqTDJohnson(FT)
    "Respiration temperature dependency for C3 plants"
    TD_R_C3::AbstractTemperatureDependency{FT} = RespirationTDCLMC3(FT)
    "Respiration temperature dependency for C4 plants"
    TD_R_C4::AbstractTemperatureDependency{FT} = RespirationTDCLMC4(FT)
    "Vcmax temperature dependency for C3 plants"
    TD_VCMAX_C3::AbstractTemperatureDependency{FT} = VcmaxTDCLMC3(FT)
    "Vcmax temperature dependency for C4 plants"
    TD_VCMAX_C4::AbstractTemperatureDependency{FT} = VcmaxTDCLMC4(FT)
    "Vpmax temperature dependency"
    TD_VPMAX::AbstractTemperatureDependency{FT} = VpmaxTDBoyd(FT)
    "Γ* temperature dependency"
    TD_Γ::AbstractTemperatureDependency{FT} = ΓStarTDCLM(FT)
    "η_C temperature dependency"
    TD_ηC::AbstractTemperatureDependency{FT} = ηCTDWang(FT)
    "η_L temperature dependency"
    TD_ηL::AbstractTemperatureDependency{FT} = ηLTDWang(FT)

    # soil albedo method
    "Soil albedo method"
    SOIL_ALBEDO::AbstractSoilAlbedo = SoilAlbedoHyperspectralCLIMA()

    # stomatal conductance model
    "Stomatal conductance model"
    STOMATAL_MODEL::AbstractStomatalConductanceModel{FT} = WangSM{FT}()
end;
