# This file contains the van der Tol et al. (2013) fluorescence model

#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2022-Jan-14: add van der Tol model struct
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
Base.@kwdef mutable struct VanDerTolFluorescenceModel{FT<:AbstractFloat}
    # General model information
    "Fitting parameter K_0"
    K_0::FT = 5.01
    "Fitting parameter α"
    K_A::FT = 1.93
    "Fitting parameter β"
    K_B::FT = 10
end

""" VanDerTolFluorescenceModel that uses data from all observations """
VDTModelAll(FT) = VanDerTolFluorescenceModel{FT}(K_0 = 2.48, K_A = 2.83, K_B = 0.114)

""" VanDerTolFluorescenceModel that uses data from drought stressed observations """
VDTModelDrought(FT) = VanDerTolFluorescenceModel{FT}(K_0 = 5.01, K_A = 1.93, K_B = 10);
