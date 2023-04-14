#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2021-Agu-10: add constructors within the structure rather than initialize it externally
#     2021-Sep-30: rename LeafBio to LeafBiophysics to be more specific
#     2021-Oct-19: sort variable to prognostic and dignostic catergories
#     2021-Oct-21: rename f_sense and K_SENES to brown and K_BROWN
#     2021-Nov-24: tease apart the characteristic absorption curves to HyperspectralAbsorption
#     2022-Jul-20: use kwdef for the constructor
#     2022-Jul-20: add field DATASET to struct
#     2022-Apr-14: add DIM_WL to HyperspectralAbsorption type
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Immutable struct that contains leaf biophysical traits used to run leaf reflection and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct HyperspectralAbsorption{FT<:AbstractFloat,DIMS}
    # Constant features
    "Specific absorption coefficients of anthocynanin `[-]`"
    K_ANT::Vector{FT}
    "Specific absorption coefficients of senescent material (brown pigments) `[-]`"
    K_BROWN::Vector{FT}
    "Specific absorption coefficients of chlorophyll a and b `[-]`"
    K_CAB::Vector{FT}
    "Specific absorption coefficients of violaxanthin carotenoid `[-]`"
    K_CAR_V::Vector{FT}
    "Specific absorption coefficients of zeaxanthin carotenoid `[-]`"
    K_CAR_Z::Vector{FT}
    "Specific absorption coefficients of carbon-based constituents `[-]`"
    K_CBC::Vector{FT}
    "Specific absorption coefficients of water `[-]`"
    K_H₂O::Vector{FT}
    "Specific absorption coefficients of dry matter `[-]`"
    K_LMA::Vector{FT}
    "Specific absorption coefficients of protein `[-]`"
    K_PRO::Vector{FT}
    "Specific absorption coefficients of PS I and II `[-]`"
    K_PS::Vector{FT}
    "Refractive index `[-]`"
    NR::Vector{FT}
end

HyperspectralAbsorption{FT,DIMS}(dset::String) where {FT,DIMS} = (
    return HyperspectralAbsorption{FT,DIMS}(
                K_ANT   = read_nc(dset, "K_ANT"),
                K_BROWN = read_nc(dset, "K_BROWN"),
                K_CAB   = read_nc(dset, "K_CAB"),
                K_CAR_V = read_nc(dset, "K_CAR_V"),
                K_CAR_Z = read_nc(dset, "K_CAR_Z"),
                K_CBC   = read_nc(dset, "K_CBC"),
                K_H₂O   = read_nc(dset, "K_H₂O"),
                K_LMA   = read_nc(dset, "K_LMA"),
                K_PRO   = read_nc(dset, "K_PRO"),
                K_PS    = read_nc(dset, "K_PS"),
                NR      = read_nc(dset, "NR")
    )
);
