"""
Method configuration for the SPAC model
"""
Base.@kwdef mutable struct SPACMethods{FT<:AbstractFloat}
    # Photosynthesis Methods
    "Photosynthesis method"
    PS_METHODS::PhotosynthesisMethods{FT} = PhotosynthesisMethods{FT}()

    # Fluorescence Methods
    "Fluorescence spectra method"
    FLUORESCENCE_SPECTRA_METHOD::Union{DualspectFluorescenceSpectra,FluspectFluorescenceSpectra,PlatespectFluorescenceSpectra} = PlatespectFluorescenceSpectra()

    # canopy radiative transfer method
    "Canopy radiative transfer method"
    CANOPY_RT_METHOD::UnionCanopyRTMethod = CanopyRTEmerald()
    "Soil albedo method"
    SOIL_ALBEDO::UnionSoilAlbedo = SoilAlbedoBroadbandCLIMA()

    # stomatal conductance model
    "Stomatal conductance model"
    STOMATAL_MODEL::UnionStomatalConductanceModel{FT} = WangSM{FT}()
end;
