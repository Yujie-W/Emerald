#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2024-Aug-16: add SIFMatrixFluspectMethod (N_DUB)
#     2024-Aug-16: add SIFMatrixPlatespectMethod
#     2024-Feb-09: add SIFMatrixDualspectMethod (N_DUB)
#
#######################################################################################################################################################################################################
""" Method to compute SIF matrices using the doubling method """
Base.@kwdef struct SIFMatrixDualspectMethod
    N::Int = 10
end;


""" Method to compute SIF matrices using the doubling method """
Base.@kwdef struct SIFMatrixFluspectMethod
    N::Int = 10
end;


""" Method to compute SIF matrices using the excitation-emission method """
struct SIFMatrixPlatespectMethod end;
