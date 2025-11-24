"""

$(TYPEDEF)

Consists a set of method selections for SPAC simulations

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SPACMethodSelections{FT}
    # Plant physiology methods
    "Stomatal model"
    STOMATAL_MODEL::AbstractStomataModel = WangSM{FT}()

    # Soil methods
    "Soil albedo method"
    SOIL_ALBEDO::AbstractSoilAlbedo = SoilAlbedoHyperspectralCLIMA()
end;
