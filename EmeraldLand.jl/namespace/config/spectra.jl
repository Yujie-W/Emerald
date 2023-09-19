#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2021-Nov-24: tease apart the characteristic absorption curves to HyperspectralAbsorption
#     2022-Jul-20: add field DATASET to struct
#     2023-Sep-13: add fields Φ_PSI and Φ_PSII
#     2023-Sep-19: add fields MAT_SOIL and SOLAR_RAD
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

    # Variable features for the soil
    "A matrix of characteristic curves"
    MAT_SOIL::Matrix{FT} = FT[read_nc(DATASET, "GSV_1") read_nc(DATASET, "GSV_2") read_nc(DATASET, "GSV_3") read_nc(DATASET, "GSV_4")]

    # Variable features for the radiation
    "Downwelling shortwave radiation reference spectrum"
    SOLAR_RAD::Matrix{FT} = [read_nc(DATASET, "E_DIFF") read_nc(DATASET, "E_DIR")]
end
