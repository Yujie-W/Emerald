"""

$(TYPEDEF)

Consists a set of method selections for SPAC simulations

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SPACMethodSelections{FT}
    # photosynthesis methods
    "C3 Ac method"
    C3_AC_METHOD::AbstractAcMethod = AcMethodC3VcmaxPi();
    "C3 Aj method"
    C3_AJ_METHOD::AbstractAjMethod = AjMethodC3JmaxPi();
    "C3 Ap method"
    C3_AP_METHOD::AbstractApMethod = ApMethodC3Vcmax();
    "C4 Ac method"
    C4_AC_METHOD::AbstractAcMethod = AcMethodC4Vcmax();
    "C4 Aj method"
    C4_AJ_METHOD::AbstractAjMethod = AjMethodC4JPSII();
    "C4 Ap method"
    C4_AP_METHOD::AbstractApMethod = ApMethodC4VcmaxPi();
    "Fluorescence method"
    FLUORESCENCE_METHOD::AbstractFluorescenceMethod{FT} = KNFluorescenceModel{FT}()

    # stomatal method
    "Stomatal model"
    STOMATAL_MODEL::AbstractStomataModel{FT} = WangSM{FT}()

    # soil albedo method
    "Soil albedo method"
    SOIL_ALBEDO::AbstractSoilAlbedo = SoilAlbedoHyperspectralCLIMA()
end;
