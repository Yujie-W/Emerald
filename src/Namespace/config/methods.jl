"""
Method configuration for the SPAC model
"""
Base.@kwdef mutable struct SPACMethods{FT<:AbstractFloat}
    # fluorescence methods
    "Fluorescence method for C3 plants"
    FLUORESCENCE_METHOD_C3::Union{CytochromeFluorescenceModel,KNFluorescenceModel{FT},QLFluorescenceModel{FT},QLFluorescenceModelHan{FT}} = KNFluorescenceModel{FT}()
    "Fluorescence method for C4 plants"
    FLUORESCENCE_METHOD_C4::Union{KNFluorescenceModel{FT},QLFluorescenceModel{FT},QLFluorescenceModelHan{FT}} = KNFluorescenceModel{FT}()
    "Fluorescence spectra method"
    FLUORESCENCE_SPECTRA_METHOD::Union{DualspectFluorescenceSpectra,FluspectFluorescenceSpectra,PlatespectFluorescenceSpectra} = PlatespectFluorescenceSpectra()

    # photosynthesis models - Ac, Aj, Ap methods
    "C3 Ac method"
    C3_AC_METHOD::AcMethodC3VcmaxPi = AcMethodC3VcmaxPi()
    "C3 Aj method"
    C3_AJ_METHOD::Union{AjMethodC3JmaxPi,AjMethodC3VqmaxPi} = AjMethodC3JmaxPi()
    "C3 Ap method"
    C3_AP_METHOD::Union{ApMethodC3Inf,ApMethodC3Vcmax} = ApMethodC3Vcmax()
    "C4 Ac method"
    C4_AC_METHOD::AcMethodC4Vcmax = AcMethodC4Vcmax()
    "C4 Aj method"
    C4_AJ_METHOD::AjMethodC4JPSII = AjMethodC4JPSII()
    "C4 Ap method"
    C4_AP_METHOD::Union{ApMethodC4VcmaxPi,ApMethodC4VpmaxPi} = ApMethodC4VcmaxPi()

    # photosynthesis models - colimitation methods
    "Colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::Union{MinimumColimit,QuadraticColimit{FT},SerialColimit,SquareColimit} = MinimumColimit()
    "Colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::Union{MinimumColimit,QuadraticColimit{FT},SerialColimit,SquareColimit} = MinimumColimit()
    "Colimitation method for J (for C3 only)"
    COLIMIT_J::Union{MinimumColimit,QuadraticColimit{FT},SerialColimit,SquareColimit} = ColimitJCLM(FT)

    # photosynthesis models - temperature dependency methods
    "Jmax temperature dependency"
    TD_JMAX::UnionTemperatureDependency{FT} = JmaxTDCLM(FT)
    "Kc temperature dependency"
    TD_KC::UnionTemperatureDependency{FT} = KcTDCLM(FT)
    "Ko temperature dependency"
    TD_KO::UnionTemperatureDependency{FT} = KoTDCLM(FT)
    "Kpep temperature dependency"
    TD_KPEP::UnionTemperatureDependency{FT} = KpepTDBoyd(FT)
    "Kpep temperature dependency to use with C4CLM method"
    TD_KPEP_CLM::UnionTemperatureDependency{FT} = Q10TDKpepCLM(FT)
    "Kq temperature dependency"
    TD_KQ::UnionTemperatureDependency{FT} = KqTDJohnson(FT)
    "Respiration temperature dependency for C3 plants"
    TD_R_C3::UnionTemperatureDependency{FT} = RespirationTDCLMC3(FT)
    "Respiration temperature dependency for C4 plants"
    TD_R_C4::UnionTemperatureDependency{FT} = RespirationTDCLMC4(FT)
    "Vcmax temperature dependency for C3 plants"
    TD_VCMAX_C3::UnionTemperatureDependency{FT} = VcmaxTDCLMC3(FT)
    "Vcmax temperature dependency for C4 plants"
    TD_VCMAX_C4::UnionTemperatureDependency{FT} = VcmaxTDCLMC4(FT)
    "Vpmax temperature dependency"
    TD_VPMAX::UnionTemperatureDependency{FT} = VpmaxTDBoyd(FT)
    "Γ* temperature dependency"
    TD_Γ::UnionTemperatureDependency{FT} = ΓStarTDCLM(FT)
    "η_C temperature dependency"
    TD_ηC::UnionTemperatureDependency{FT} = ηCTDWang(FT)
    "η_L temperature dependency"
    TD_ηL::UnionTemperatureDependency{FT} = ηLTDWang(FT)

    # canopy radiative transfer method
    "Canopy radiative transfer method"
    CANOPY_RT_METHOD::UnionCanopyRTMethod = CanopyRTEmerald()
    "Soil albedo method"
    SOIL_ALBEDO::UnionSoilAlbedo = SoilAlbedoHyperspectralCLIMA()

    # stomatal conductance model
    "Stomatal conductance model"
    STOMATAL_MODEL::UnionStomatalConductanceModel{FT} = WangSM{FT}()
end;
