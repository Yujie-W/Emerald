# The file contains the definitions of the pressure volume curve models
# Because PV curve is a trait (a state), so there is no need to split the struct into state and auxilary parts.

#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Apr-20: add abstract type for pressure volume curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractPVCurve:
- [`LinearPVCurve`](@ref)
- [`SegmentedPVCurve`](@ref)

"""
abstract type AbstractPVCurve{FT<:AbstractFloat} end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-22: add exponential PV curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains information for exponential PV curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ExponentialPVCurve{FT<:AbstractFloat} <: AbstractPVCurve{FT}
    # General model information
    "Conductance for refilling (relative to maximum) `[MPa⁻¹ s⁻¹]`"
    k_refill::FT = 1e-4
    "Slope of the linear PV curve (relative to maximum) `[MPa⁻¹]`"
    slope::FT = 0.2
end;

sync_state!(state_from::ExponentialPVCurve{FT}, state_to::ExponentialPVCurve{FT}) where {FT} = (
    state_to.k_refill = state_from.k_refill;
    state_to.slope = state_from.slope;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Apr-20: add linear PV curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains information for linear PV curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LinearPVCurve{FT<:AbstractFloat} <: AbstractPVCurve{FT}
    # General model information
    "Conductance for refilling (relative to maximum) `[MPa⁻¹ s⁻¹]`"
    k_refill::FT = 1e-4
    "Slope of the linear PV curve (relative to maximum) `[MPa⁻¹]`"
    slope::FT = 0.2
end;

sync_state!(state_from::LinearPVCurve{FT}, state_to::LinearPVCurve{FT}) where {FT} = (
    state_to.k_refill = state_from.k_refill;
    state_to.slope = state_from.slope;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-24: add segmented PV curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains information for segmented PV curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SegmentedPVCurve{FT<:AbstractFloat} <: AbstractPVCurve{FT}
    # General model information
    "n_o / maximum V `[mol m⁻³]`"
    c_all::FT = 300
    "Conductance for refilling (relative to maximum) `[MPa⁻¹ s⁻¹]`"
    k_refill::FT = 1e-4
    "Apoplastic water content relative to maximum water volume"
    θ_apo::FT = 0.4
    "Relative water content at turgor loss point"
    θ_tlp::FT = 0.8
    "Bulk modulus of elasticity `[MPa]`"
    ϵ_bulk::FT = 20
end;

sync_state!(state_from::SegmentedPVCurve{FT}, state_to::SegmentedPVCurve{FT}) where {FT} = (
    state_to.c_all = state_from.c_all;
    state_to.k_refill = state_from.k_refill;
    state_to.θ_apo = state_from.θ_apo;
    state_to.θ_tlp = state_from.θ_tlp;
    state_to.ϵ_bulk = state_from.ϵ_bulk;

    return nothing
);
