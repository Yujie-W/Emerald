# This file contains the spectra information to use with the SPAC system
# By default, the reference spectra is to run at hyperspectral mode
# However, if is also possible to run the model at broadband mode, and we provide a constructor for this purpose (SIF will be disabled this way)

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2022-Jul-20: add field DATASET to struct
#     2023-Sep-11: add field ΔΛ_SIF
#     2023-Sep-13: add fields Φ_PSI and Φ_PSII
#     2023-Sep-19: add fields MAT_SOIL and SOLAR_RAD
#     2023-Oct-17: add field ρ_STEM
#     2024-Jan-23: change default minimum PAR to 300 nm
#     2024-Jan-23: add field for PAR using the 300-700 nm definition
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Immutable struct that contains leaf biophysical traits used to run leaf reflection and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct ReferenceSpectra{FT<:AbstractFloat}
    # File path to the Netcdf dataset
    "File path to the Netcdf dataset"
    DATASET::String = LAND_2021

    # Wavelegnth bins and upper/lower boundaries
    "Wavelength (bins) `[nm]`"
    Λ::Vector{FT} = read_nc(DATASET, "WL")
    "Lower boundary wavelength `[nm]`"
    Λ_LOWER::Vector{FT} = read_nc(DATASET, "WL_LOWER")
    "Upper boundary wavelength `[nm]`"
    Λ_UPPER::Vector{FT} = read_nc(DATASET, "WL_UPPER")
    "Differential wavelength `[nm]`"
    ΔΛ::Vector{FT} = Λ_UPPER .- Λ_LOWER

    # Constant features for the leaf
    "Specific absorption coefficients of anthocynanin `[-]`"
    K_ANT::Vector{FT} = read_nc(DATASET, "K_ANT")
    "Specific absorption coefficients of senescent material (brown pigments) `[-]`"
    K_BROWN::Vector{FT} = read_nc(DATASET, "K_BROWN")
    "Specific absorption coefficients of chlorophyll a and b `[-]`"
    K_CAB::Vector{FT} = read_nc(DATASET, "K_CAB")
    "Specific absorption coefficients of violaxanthin carotenoid `[-]`"
    K_CAR_V::Vector{FT} = read_nc(DATASET, "K_CAR_V")
    "Specific absorption coefficients of zeaxanthin carotenoid `[-]`"
    K_CAR_Z::Vector{FT} = read_nc(DATASET, "K_CAR_Z")
    "Specific absorption coefficients of carbon-based constituents `[-]`"
    K_CBC::Vector{FT} = read_nc(DATASET, "K_CBC")
    "Specific absorption coefficients of water `[-]`"
    K_H₂O::Vector{FT} = read_nc(DATASET, "K_H₂O")
    "Specific absorption coefficients of dry matter `[-]`"
    K_LMA::Vector{FT} = read_nc(DATASET, "K_LMA")
    "Specific absorption coefficients of protein `[-]`"
    K_PRO::Vector{FT} = read_nc(DATASET, "K_PRO")
    "Refractive index `[-]`"
    NR::Vector{FT} = read_nc(DATASET, "NR")
    "Fluorescence yield of PS I and II probability function `[nm⁻¹]`"
    Φ_PS::Vector{FT} = read_nc(DATASET, "K_PS")
    "Fluorescence yield of PS I probability function `[nm⁻¹]`"
    Φ_PSI::Vector{FT} = read_nc(DATASET, "K_PS1")
    "Fluorescence yield of PS II probability function `[nm⁻¹]`"
    Φ_PSII::Vector{FT} = read_nc(DATASET, "K_PS2")

    # Stem reflectance
    "Stem reflectance `[-]`"
    ρ_STEM::Vector{FT} = zeros(FT, length(Λ)) .+ 0.2

    # Variable features for the soil
    "A matrix of characteristic curves"
    MAT_SOIL::Matrix{FT} = FT[read_nc(DATASET, "GSV_1") read_nc(DATASET, "GSV_2") read_nc(DATASET, "GSV_3") read_nc(DATASET, "GSV_4")]

    # Variable features for the radiation
    "Downwelling shortwave radiation reference spectrum"
    SOLAR_RAD::Matrix{FT} = [read_nc(DATASET, "E_DIR") read_nc(DATASET, "E_DIFF")]

    # Constants
    "Wavelength limits for NIR `[nm]`"
    WL_NIR::Vector{FT} = FT[700, 2500]
    "Wavelength limits for PAR `[nm]`"
    WL_PAR::Vector{FT} = FT[300, 750]
    "Wavelength limits for PAR using the 300-700 nm setup `[nm]`"
    WL_PAR_700::Vector{FT} = FT[300, 700]
    "Wavelength limits for SIF emission `[nm]`"
    WL_SIF::Vector{FT} = FT[640, 850]
    "Wavelength limits for SIF excitation `[nm]`"
    WL_SIFE::Vector{FT} = FT[300, 750]

    # Indices
    "Indicies of Λ_NIR in Λ"
    IΛ_NIR::Vector{Int} = findall( WL_NIR[1] .<= Λ .<= WL_NIR[2] )
    "Indicies of Λ_PAR in Λ"
    IΛ_PAR::Vector{Int} = findall( WL_PAR[1] .<= Λ .<= WL_PAR[2] )
    "Indicies of Λ_PAR in Λ using the 300-700 nm definition"
    IΛ_PAR_700::Vector{Int} = findall( WL_PAR_700[1] .<= Λ .<= WL_PAR_700[2] )
    "Indicies of Λ_SIF in Λ"
    IΛ_SIF::Vector{Int} = findall( WL_SIF[1] .<= Λ .<= WL_SIF[2] )
    "Indicies of Λ_SIFE in Λ"
    IΛ_SIFE::Vector{Int} = findall( WL_SIFE[1] .<= Λ .<= WL_SIFE[2] )

    # Constants based on the ones above
    "Differential wavelength for PAR `[nm]`"
    ΔΛ_PAR::Vector{FT} = ΔΛ[IΛ_PAR]
    "Differential wavelength for PAR using 300-700 nm definition `[nm]`"
    ΔΛ_PAR_700::Vector{FT} = ΔΛ[IΛ_PAR_700]
    "Differential wavelength for SIF `[nm]`"
    ΔΛ_SIF::Vector{FT} = ΔΛ[IΛ_SIF]
    "Differential wavelength for SIF excitation `[nm]`"
    ΔΛ_SIFE::Vector{FT} = ΔΛ[IΛ_SIFE]
    "Wavelength bins for PAR `[nm]`"
    Λ_PAR::Vector{FT} = Λ[IΛ_PAR]
    "Wavelength bins for PAR using the 300-700 nm definition `[nm]`"
    Λ_PAR_700::Vector{FT} = Λ[IΛ_PAR_700]
    "Wavelength bins for SIF `[nm]`"
    Λ_SIF::Vector{FT} = Λ[IΛ_SIF]
    "Wavelength bins for SIF excitation `[nm]`"
    Λ_SIFE::Vector{FT} = Λ[IΛ_SIFE]
end;


"""

    broadband_spectra(FT)

Return a ReferenceSpectra object with broadband spectra, given
- `FT`: the floating point type to use

"""
broadband_spectra(FT) = ReferenceSpectra{FT}(
            DATASET   = "Customized Broadband Spectra",
            Λ         = [550, 1600],
            Λ_LOWER   = [400, 700],
            Λ_UPPER   = [700, 2500],
            ΔΛ        = [300, 1800],
            K_ANT     = [0, 0],                             # will not be used
            K_BROWN   = [0, 0],                             # will not be used
            K_CAB     = [0, 0],                             # will not be used
            K_CAR_V   = [0, 0],                             # will not be used
            K_CAR_Z   = [0, 0],                             # will not be used
            K_CBC     = [0, 0],                             # will not be used
            K_H₂O     = [0, 0],                             # will not be used
            K_LMA     = [0, 0],                             # will not be used
            K_PRO     = [0, 0],                             # will not be used
            NR        = [1.5, 1.5],                         # will not be used
            Φ_PS      = [0, 0],                             # will not be used
            Φ_PSI     = [0, 0],                             # will not be used
            Φ_PSII    = [0, 0],                             # will not be used
            ρ_STEM    = [0.16, 0.39],
            MAT_SOIL  = [0 0; 0 0],                         # will not be used
            SOLAR_RAD = [191.544 119.923; 277.752 55.076],
            WL_NIR    = [700, 2500],
            WL_PAR    = [400, 700],
            WL_SIF    = [640, 850],                         # will not be used
            WL_SIFE   = [400, 750],                         # will not be used
            IΛ_NIR    = [2],
            IΛ_PAR    = [1],
            IΛ_SIF    = [1],                                # will not be used
            IΛ_SIFE   = [1],                                # will not be used
            ΔΛ_PAR    = [300],
            ΔΛ_SIF    = [210],                              # will not be used
            ΔΛ_SIFE   = [350],                              # will not be used
            Λ_PAR     = [550],
            Λ_SIF     = [745],                              # will not be used
            Λ_SIFE    = [550],                              # will not be used
);
