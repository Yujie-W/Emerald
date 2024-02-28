# This file contains the capacitor outside the xylem, e.g., leaf capacitor and root-trunk junction capacitor

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add struct ExtraXylemCapacitorTrait
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for extraxylary capacitor

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ExtraXylemCapacitorTrait{FT}
    "Maximum extraxylary hydraulic conductance `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_max::FT = 100
    "Pressure volume curve of the capacitor"
    pv::Union{ExponentialPVCurve{FT}, LinearPVCurve{FT}, SegmentedPVCurve{FT}} = SegmentedPVCurve{FT}()
    "Capacitor maximum volume per basal area or per leaf area `[mol m⁻²]`"
    v_max::FT = 5
    "Vulnerability curve of the extraxylary capacitor"
    vc::Union{ComplexVC{FT}, LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}} = WeibullVC{FT}(5,1)
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: define the struct to store the state variables used in extraxylary capacitor
#     2024-Feb-28: move field p_leaf to state (see the comment in the state struct for the reason)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables for extraxylary capacitor

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ExtraXylemCapacitorState{FT}
    "Current volume of the capacitor `[mol]`"
    v_storage::FT = 0
    # note: This in theory should not be a state variable, but it is used to compute saturation vapor pressure in the next time step and then flow profile.
    #       However, if we store this as an auxiliary variable, then we need to compute it after we know the flow profile. It is a chicken and egg problem, so we store it as a state here.
    #       But, we still compute it still as an auxiliary variable.
    "Pressure at the end; of the capacitor `[MPa]`"
    p_leaf::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: define the struct to store the auxiliary variables used in extraxylary capacitor
#     2023-Sep-27: add field p
#     2024-Feb-28: move field p_leaf to state (see the comment in the state struct for the reason)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxiliary variables for extraxylary capacitor

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ExtraXylemCapacitorAuxil{FT}
    "Flow rate out from the capacitor `[mol s⁻¹]`"
    flow::FT = 0
    "Hydraulic conductance of the capacitor `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k::FT = 100
    "pressure of the capacitor `[MPa]`"
    p::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: define the struct to store the capacitor
#     2024-Feb-26: add field trait
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the capacitor state and auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ExtraXylemCapacitor{FT}
    "Trait variables of the capacitor"
    trait::ExtraXylemCapacitorTrait{FT} = ExtraXylemCapacitorTrait{FT}()
    "State variables of the capacitor"
    state::ExtraXylemCapacitorState{FT} = ExtraXylemCapacitorState{FT}()
    "Auxilary variables of the capacitor"
    auxil::ExtraXylemCapacitorAuxil{FT} = ExtraXylemCapacitorAuxil{FT}()
end;
