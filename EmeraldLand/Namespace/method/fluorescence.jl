# This file contains the van der Tol et al. (2013) fluorescence model

#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2022-Jan-14: add Kn based model
# Sources
#     van der Tol et al. (2013) Models of fluorescence and photosynthesis for interpreting measurements of solar-induced chlorophyll fluorescence
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores van der Tol et al. (2013) fluorescence model parameters.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct KNFluoscenceModel{FT<:AbstractFloat}
    # General model information
    "Fitting parameter K_0"
    K_0::FT = 5.01
    "Fitting parameter α"
    K_A::FT = 1.93
    "Fitting parameter β"
    K_B::FT = 10
end;

KNFluoscenceModelAll(FT) = KNFluoscenceModel{FT}(K_0 = 2.48, K_A = 2.83, K_B = 0.114)

KNFluoscenceModelDrought(FT) = KNFluoscenceModel{FT}(K_0 = 5.01, K_A = 1.93, K_B = 10);


#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2023-Oct-27: add qL based model
#     2023-Oct-30: remove K_A from the method (which should be 1)
# Sources
#     Han et al. (2022) The physiological basis for estimating photosynthesis from Chla fluorescence
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores Han et al. (2022) fluorescence model parameters.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct QLFluoscenceModel{FT<:AbstractFloat}
    "Fitting parameter qb"
    K_B::FT = 0.95e-3
end;

QLFluoscenceModelC3(FT) = QLFluoscenceModel{FT}(K_B = 0.95e-3);

QLFluoscenceModelC4(FT) = QLFluoscenceModel{FT}(K_B = 0.63e-3);