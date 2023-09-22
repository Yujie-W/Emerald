# This file contains the state and auxilary variables for xylem hydraulics

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-22: define the struct to store the state variables used in xylem hydraulics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables for xylem hydraulics

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct XylemHydraulicsState{FT}
    "Area of xylem (root and stem) or of leaf `[m²]`"
    area::FT = 1
    "Maximal xylem hydraulic conductivity `[mol s⁻¹ MPa⁻¹ m⁻¹]` for root and stem; `[mol s⁻¹ MPa⁻¹ m⁻²]` for leaf"
    k_max::FT = 25
    "Length `[m]`"
    l::FT = 1
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT}
    "Storage per element `[mol]`"
    v_storage::Vector{FT}
    "Vulnerability curve"
    vc::AbstractXylemVC{FT} = WeibullVC{FT}()
    "Height change `[m]`"
    Δh::FT = 1
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-22: define the struct to store the non-steady state auxilary variables used in xylem hydraulics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxilary variables for xylem hydraulics at non-steady state mode

# Fields

$(TYPEDFIELDS)

"""
mutable struct XylemHydraulicsAuxilNSS{FT} <: AbstractFlowProfile{FT}
    "Critical flow rate `[mol s⁻¹ m⁻²]`"
    e_crit::FT
    "Flow rates at each plane (N+1) `[mol s⁻¹]` for root and stem or `[mol m⁻² s⁻¹]` for leaf"
    flow::Vector{FT}
    "Flow rates from the buffer system at each slice (N) `[mol s⁻¹]` for root and stem or `[mol m⁻² s⁻¹]` for leaf"
    flow_buffer::Vector{FT}
    "Vector of leaf kr history per element"
    k_history::Vector{FT}
    "Pressure of storage per element"
    p_storage::Vector{FT}
    "Vector of xylem water pressure at each plance (N+1) `[MPa]`"
    ps::Vector{FT}
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-22: define the struct to store the steady state auxilary variables used in xylem hydraulics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxilary variables for xylem hydraulics at steady state mode

# Fields

$(TYPEDFIELDS)

"""
mutable struct XylemHydraulicsAuxilSS{FT} <: AbstractFlowProfile{FT}
    "Critical flow rate `[mol s⁻¹ m⁻²]`"
    e_crit::FT
    "Flow rates through the system `[mol s⁻¹]` for root and stem or `[mol m⁻² s⁻¹]` for leaf"
    flow::FT
    "Vector of leaf kr history per element"
    k_history::Vector{FT}
    "Vector of xylem water pressure at each plance (N+1) `[MPa]`"
    ps::Vector{FT}
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-22: define the struct to store the auxilary variables used in xylem hydraulics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state and auxilary fields for xylem hydraulics

# Fields

$(TYPEDFIELDS)

"""
mutable struct XylemHydraulics{FT}
    "State variables"
    state::XylemHydraulicsState{FT}
    "Auxiliary variables"
    auxil::Union{XylemHydraulicsAuxilNSS{FT}, XylemHydraulicsAuxilSS{FT}}
end;
