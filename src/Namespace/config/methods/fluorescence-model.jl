"""
Module that defines fluorescence methods used in leaf photosynthesis modeling:
"""
abstract type AbstractFluorescenceMethod{FT<:AbstractFloat} end;


"""
Structure for the C3 Cytochrome fluorescence model.
"""
struct CytochromeFluorescenceModel{FT<:AbstractFloat} <: AbstractFluorescenceMethod{FT} end;


"""
Structure that stores van der Tol et al. (2013) fluorescence model parameters.
- van der Tol et al. (2013) Models of fluorescence and photosynthesis for interpreting measurements of solar-induced chlorophyll fluorescence
"""
Base.@kwdef mutable struct KNFluorescenceModel{FT<:AbstractFloat} <: AbstractFluorescenceMethod{FT}
    # General model information
    "Fitting parameter K_0"
    K_0::FT = 5.01
    "Fitting parameter α"
    K_A::FT = 1.93
    "Fitting parameter β"
    K_B::FT = 10
end;


"""
Structure that stores modified Han et al. (2022) fluorescence model parameters.
- Han et al. (2022) The physiological basis for estimating photosynthesis from Chla fluorescence
"""
Base.@kwdef mutable struct QLFluorescenceModel{FT<:AbstractFloat} <: AbstractFluorescenceMethod{FT}
    "Fitting parameter qb"
    K_B::FT = 0.95e-3 / 0.85
end;


"""
Structure that stores original Han et al. (2022) fluorescence model parameters.
- Han et al. (2022) The physiological basis for estimating photosynthesis from Chla fluorescence
"""
Base.@kwdef mutable struct QLFluorescenceModelHan{FT<:AbstractFloat} <: AbstractFluorescenceMethod{FT}
    "Fitting parameter α"
    K_A::FT = 0.8
    "Fitting parameter β"
    K_B::FT = 0.95e-3 / 0.85
end;
