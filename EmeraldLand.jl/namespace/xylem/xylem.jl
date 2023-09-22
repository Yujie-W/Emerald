# This file contains the state and auxilary variables for xylem hydraulics

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-22: define the struct to store the state variables used in xylem hydraulics
#     2023-Sep-22: add field v_max, pv
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
    "Pressure volume curve"
    pv::AbstractPVCurve{FT} = LinearPVCurve{FT}()
    "Storage per element `[mol]`"
    v_storage::Vector{FT}
    "Maximum capaciatance per volume of wood `[mol m⁻³]`"
    v_max::FT
    "Vulnerability curve"
    vc::AbstractXylemVC{FT} = WeibullVC{FT}()
    "Height change `[m]`"
    Δh::FT = 1
end;

XylemHydraulicsState(config::SPACConfiguration{FT}; area::Number = 1, l::Number = 1, v_max::Number = 0.1 * ρ_H₂O() / M_H₂O()) where {FT} = XylemHydraulicsState{FT}(
            area      = area,
            l         = l,
            p_history = zeros(FT, config.DIM_XYLEM),
            v_storage = ones(FT, config.DIM_XYLEM) ./ config.DIM_XYLEM .* (v_max * area * l),
            v_max     = v_max
);


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
Base.@kwdef mutable struct XylemHydraulicsAuxilNSS{FT} <: AbstractFlowProfile{FT}
    "Critical flow rate `[mol s⁻¹ m⁻²]`"
    e_crit::FT = 0
    "Flow rates at each plane (N+1) `[mol s⁻¹]` for root and stem or `[mol m⁻² s⁻¹]` for leaf"
    flow::Vector{FT}
    "Flow rates from the buffer system at each slice (N) `[mol s⁻¹]` for root and stem or `[mol m⁻² s⁻¹]` for leaf"
    flow_buffer::Vector{FT}
    "Vector of leaf kr history per element"
    k_history::Vector{FT}
    "Pressure of storage per element"
    p_storage::Vector{FT}
    "Vector of xylem water pressure at each plance (N+1) `[MPa]`"
    pressure::Vector{FT}
end;

XylemHydraulicsAuxilNSS(config::SPACConfiguration{FT}) where {FT} = XylemHydraulicsAuxilNSS{FT}(
            flow        = zeros(FT, config.DIM_XYLEM + 1),
            flow_buffer = zeros(FT, config.DIM_XYLEM),
            k_history   = zeros(FT, config.DIM_XYLEM),
            p_storage   = zeros(FT, config.DIM_XYLEM),
            pressure    = zeros(FT, config.DIM_XYLEM + 1)
);


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
Base.@kwdef mutable struct XylemHydraulicsAuxilSS{FT} <: AbstractFlowProfile{FT}
    "Critical flow rate `[mol s⁻¹ m⁻²]`"
    e_crit::FT = 0
    "Flow rates through the system `[mol s⁻¹]` for root and stem or `[mol m⁻² s⁻¹]` for leaf"
    flow::FT = 0
    "Vector of leaf kr history per element"
    k_history::Vector{FT}
    "Vector of xylem water pressure at each plance (N+1) `[MPa]`"
    pressure::Vector{FT}
end;

XylemHydraulicsAuxilSS(config::SPACConfiguration{FT}) where {FT} = XylemHydraulicsAuxilSS{FT}(
            k_history = zeros(FT, config.DIM_XYLEM),
            pressure  = zeros(FT, config.DIM_XYLEM + 1)
);


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
Base.@kwdef mutable struct XylemHydraulics{FT}
    "State variables"
    state::XylemHydraulicsState{FT}
    "Auxiliary variables"
    auxil::Union{XylemHydraulicsAuxilNSS{FT}, XylemHydraulicsAuxilSS{FT}}
end;

XylemHydraulics(config::SPACConfiguration{FT}) where {FT} = XylemHydraulics{FT}(
            state = XylemHydraulicsState(config),
            auxil = config.STEADY_STATE_FLOW ? XylemHydraulicsAuxilSS(config) : XylemHydraulicsAuxilNSS(config)
);
