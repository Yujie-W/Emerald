"""
Hierarchy of methods to compute fluorescence spectra
- DualspectFluorescenceSpectra
- FluspectFluorescenceSpectra
- PlatespectFluorescenceSpectra
"""
abstract type AbstractFluorescenceSpectraMethod end;


""" Method to compute SIF matrices using the doubling method """
Base.@kwdef struct DualspectFluorescenceSpectra <: AbstractFluorescenceSpectraMethod
    N::Int = 10
end;


""" Method to compute SIF matrices using the doubling method """
Base.@kwdef struct FluspectFluorescenceSpectra <: AbstractFluorescenceSpectraMethod
    N::Int = 10
end;


""" Method to compute SIF matrices using the excitation-emission method """
struct PlatespectFluorescenceSpectra <: AbstractFluorescenceSpectraMethod end;
