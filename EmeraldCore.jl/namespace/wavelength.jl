#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-10: refactor the structure with renamed fields
#     2021-Aug-10: add a constructor within the structure to avoid external initialization
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2022-Jul-20: use kwdef for the constructor
#     2022-Jul-20: add field DATASET to struct
#     2022-Apr-14: add DIM_WL to WaveLengthSet type
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Immutable structure that stores wave length information.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct WaveLengthSet{FT,DIMS}
    # Constants
    "Wavelength (bins) `[nm]`"
    Λ::Vector{FT}
    "Lower boundary wavelength `[nm]`"
    Λ_LOWER::Vector{FT}
    "Upper boundary wavelength `[nm]`"
    Λ_UPPER::Vector{FT}

    # Indices
    "Indicies of Λ_NIR in Λ"
    IΛ_NIR::Vector{Int}
    "Indicies of Λ_PAR in Λ"
    IΛ_PAR::Vector{Int}
    "Indicies of Λ_SIF in Λ"
    IΛ_SIF::Vector{Int}
    "Indicies of Λ_SIFE in Λ"
    IΛ_SIFE::Vector{Int}

    # Constants based on the ones above
    "Differential wavelength `[nm]`"
    ΔΛ::Vector{FT} = Λ_UPPER .- Λ_LOWER
    "Differential wavelength for PAR `[nm]`"
    ΔΛ_PAR::Vector{FT} = ΔΛ[IΛ_PAR]
    "Differential wavelength for SIF excitation `[nm]`"
    ΔΛ_SIFE::Vector{FT} = ΔΛ[IΛ_SIFE]
    "Wavelength bins for PAR `[nm]`"
    Λ_PAR::Vector{FT} = Λ[IΛ_PAR]
    "Wavelength bins for SIF `[nm]`"
    Λ_SIF::Vector{FT} = Λ[IΛ_SIF]
    "Wavelength bins for SIF excitation `[nm]`"
    Λ_SIFE::Vector{FT} = Λ[IΛ_SIFE]
end

WaveLengthSet{FT,DIMS}(gcf::GeneralConfiguration) where {FT,DIMS} =  (
    _λ = read_nc(gcf.DATASET, "WL");

    return WaveLengthSet{FT,DIMS}(
                Λ       = _λ,
                Λ_LOWER = read_nc(gcf.DATASET, "WL_LOWER"),
                Λ_UPPER = read_nc(gcf.DATASET, "WL_UPPER"),
                IΛ_NIR  = findall(gcf.WL_NIR[1]  .<= _λ .<= gcf.WL_NIR[2]),
                IΛ_PAR  = findall(gcf.WL_PAR[1]  .<= _λ .<= gcf.WL_PAR[2]),
                IΛ_SIF  = findall(gcf.WL_SIF[1]  .<= _λ .<= gcf.WL_SIF[2]),
                IΛ_SIFE = findall(gcf.WL_SIFE[1] .<= _λ .<= gcf.WL_SIFE[2]),
    )
);
