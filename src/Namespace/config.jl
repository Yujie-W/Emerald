"""
General configuration for the SPAC model.
"""
Base.@kwdef mutable struct SPACConfig{FT<:AbstractFloat}
    "General configuration information"
    CONFIG_INFO::SPACConfigInfo = SPACConfigInfo()
    "Constants used in the SPAC model"
    CONSTANTS::SPACConstants{FT} = SPACConstants{FT}()
    "Dimensions of the SPAC system"
    DIMENSIONS::SPACDimensions{FT} = SPACDimensions{FT}()
    "Features on/off/settings of the SPAC model"
    FEATURES::SPACFeatures = SPACFeatures()
    "Methods used in the SPAC model"
    METHODS::SPACMethods{FT} = SPACMethods{FT}()
end;

"""

    SPACConfig(FT::DataType;
               dataset::String = OLD_PHI_2021,
               jld2_file::String = LAND_ARTIFACT,
               wl_par::Vector = [300,750],
               wl_par_700::Vector = [300,700],
               wl_selection::Union{Nothing,Vector} = nothing)

Create and return a SPAC configuration, given
- `FT` the floating number type
- `dataset` the dataset name in the JLD2 file
- `jld2_file` the JLD2 file name
- `wl_par` the wavelength range for PAR
- `wl_par_700` the wavelength range for PAR 700
- `wl_selection` the wavelength selection

"""
SPACConfig(FT::DataType;
           dataset::String = OLD_PHI_2021,
           jld2_file::String = LAND_ARTIFACT,
           wl_par::Vector = [300,750],
           wl_par_700::Vector = [300,700],
           wl_selection::Union{Nothing,Vector} = nothing) = (
    return SPACConfig{FT}(
        CONFIG_INFO = SPACConfigInfo(DATASET = dataset, JLD2_FILE = jld2_file),
        CONSTANTS   = SPACConstants{FT}(SPECTRA = ReferenceSpectra{FT}(jld2_file, dataset; wl_par = wl_par, wl_par_700 = wl_par_700, wl_selection = wl_selection)),
    )
);
