# This file contains the state and auxiliary variables for xylem hydraulics

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add struct XylemHydraulicsTrait
#     2024-Aug-29: add xylem respiration rate
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for xylem hydraulics

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct XylemHydraulicsTrait{FT}
    "Area of xylem (root and stem) or of leaf `[m²]`"
    area::FT = 1
    "Length `[m]`"
    l::FT = 1
    "Height change `[m]`"
    Δh::FT = 1

    "Heat capacity of the xylem per volume or of leaf per mass `[J m⁻³ K⁻¹]`"
    cp::FT = 1e6
    "Wood density `[kg m⁻³]`"
    ρ::FT = 600

    "Maximal xylem hydraulic conductivity `[mol s⁻¹ MPa⁻¹ m⁻¹]` for root and stem; `[mol s⁻¹ MPa⁻¹ m⁻²]` for leaf"
    k_max::FT = 100
    "Pressure volume curve"
    pv::Union{ExponentialPVCurve{FT}, LinearPVCurve{FT}, SegmentedPVCurve{FT}} = LinearPVCurve{FT}()
    "Maximum capaciatance per volume of wood `[mol m⁻³]`"
    v_max::FT = 0.1 * ρ_H₂O() / M_H₂O()
    "Vulnerability curve"
    vc::Union{ComplexVC{FT}, LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}} = WeibullVC{FT}()

    # TODO: move this out in the future to another struct
    "Respiration rate"
    r_wood::Q10{FT} = Q10TDAngiosperm(FT)
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-22: define the struct to store the state variables used in xylem hydraulics
#     2023-Sep-22: add field v_max, pv
#     2023-Sep-23: add field cp
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables for xylem hydraulics

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct XylemHydraulicsState{FT}
    "Sap wood area"
    asap::FT = 1
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT}
    "Storage per element `[mol]`"
    v_storage::Vector{FT}
end;

XylemHydraulicsState(config::SPACConfiguration{FT}) where {FT} = XylemHydraulicsState{FT}(
            p_history = zeros(FT, config.DIM_XYLEM),
            v_storage = ones(FT, config.DIM_XYLEM) ./ config.DIM_XYLEM .* (0.1 * ρ_H₂O() / M_H₂O()),
);

kill_plant!(st::XylemHydraulicsState{FT}) where {FT} = (
    st.asap = 0;
    st.p_history .= 0;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-22: define the struct to store the non-steady state auxiliary variables used in xylem hydraulics
#     2023-Sep-30: add field connected
#     2024-Aug-05: remove k_history from auxil (computed on the fly)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxiliary variables for xylem hydraulics at non-steady state mode

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct XylemHydraulicsAuxilNSS{FT}
    "Connected to the soil (for root only)"
    connected::Bool = true
    "Critical flow rate `[mol s⁻¹]`"
    e_crit::FT = 0
    "Flow rates at each plane (N+1) `[mol s⁻¹]` for root and stem or `[mol m⁻² s⁻¹]` for leaf"
    flow::Vector{FT}
    "Flow rates from the buffer system at each slice (N) `[mol s⁻¹]` for root and stem or `[mol m⁻² s⁻¹]` for leaf"
    flow_buffer::Vector{FT}
    "Pressure of storage per element"
    p_storage::Vector{FT}
    "Vector of xylem water pressure at each plance (N+1) `[MPa]`"
    pressure::Vector{FT}
end;

XylemHydraulicsAuxilNSS(config::SPACConfiguration{FT}) where {FT} = XylemHydraulicsAuxilNSS{FT}(
            flow        = zeros(FT, config.DIM_XYLEM + 1),
            flow_buffer = zeros(FT, config.DIM_XYLEM),
            p_storage   = zeros(FT, config.DIM_XYLEM),
            pressure    = zeros(FT, config.DIM_XYLEM + 1)
);

kill_plant!(st::XylemHydraulicsAuxilNSS{FT}) where {FT} = (
    st.connected = false;
    st.e_crit = 0;
    st.flow .= 0;
    st.flow_buffer .= 0;
    st.p_storage .= 0;
    st.pressure .= 0;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-22: define the struct to store the steady state auxiliary variables used in xylem hydraulics
#     2023-Sep-30: add field connected
#     2024-Aug-05: remove k_history from auxil (computed on the fly)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxiliary variables for xylem hydraulics at steady state mode

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct XylemHydraulicsAuxilSS{FT}
    "Connected to the soil (for root only)"
    connected::Bool = true
    "Critical flow rate `[mol s⁻¹]`"
    e_crit::FT = 0
    "Flow rates through the system `[mol s⁻¹]` for root and stem or `[mol m⁻² s⁻¹]` for leaf"
    flow::FT = 0
    "Vector of xylem water pressure at each plance (N+1) `[MPa]`"
    pressure::Vector{FT}
end;

XylemHydraulicsAuxilSS(config::SPACConfiguration{FT}) where {FT} = XylemHydraulicsAuxilSS{FT}(pressure  = zeros(FT, config.DIM_XYLEM + 1));

kill_plant!(st::XylemHydraulicsAuxilSS{FT}) where {FT} = (
    st.connected = false;
    st.e_crit = 0;
    st.flow = 0;
    st.pressure .= 0;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-22: define the struct to store the auxiliary variables used in xylem hydraulics
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state and auxiliary fields for xylem hydraulics

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct XylemHydraulics{FT}
    "Trait variables"
    trait::XylemHydraulicsTrait{FT} = XylemHydraulicsTrait{FT}()
    "State variables"
    state::XylemHydraulicsState{FT}
    "Auxiliary variables"
    auxil::Union{XylemHydraulicsAuxilNSS{FT}, XylemHydraulicsAuxilSS{FT}}
end;

XylemHydraulics(config::SPACConfiguration{FT}) where {FT} = XylemHydraulics{FT}(
            state = XylemHydraulicsState(config),
            auxil = config.STEADY_STATE_FLOW ? XylemHydraulicsAuxilSS(config) : XylemHydraulicsAuxilNSS(config)
);

kill_plant!(st::XylemHydraulics{FT}) where {FT} = (
    kill_plant!(st.trait);
    kill_plant!(st.state);
    kill_plant!(st.auxil);

    return nothing
);
