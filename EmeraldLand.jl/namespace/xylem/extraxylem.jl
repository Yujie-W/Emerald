# This file contains the capacitor outside the xylem, e.g., leaf capacitor and root-trunk junction capacitor

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: define the struct to store the state variables used in extraxylary capacitor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables for extraxylary capacitor

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ExtraXylemCapacitorState{FT}
    "Maximum extraxylary hydraulic conductance `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k_max::FT = 100
    "Pressure volume curve of the capacitor"
    pv::AbstractPVCurve{FT} = ExponentialPVCurve{FT}()
    "Current volume of the capacitor `[mol m⁻²]`"
    v_storage::FT = 0
    "Capacitor maximum volume per basal area or per leaf area `[mol m⁻²]`"
    v_max::FT = 0
    "Vulnerability curve of the extraxylary capacitor"
    vc::AbstractXylemVC{FT} = WeibullVC{FT}(5,1)
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: define the struct to store the auxilary variables used in extraxylary capacitor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxilary variables for extraxylary capacitor

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ExtraXylemCapacitorAuxil{FT}
    "Flow rate out from the capacitor `[mol m⁻² s⁻¹]`"
    flow::FT = 0
    "Hydraulic conductance of the capacitor `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    k::FT = 100
    "Pressure at the end of the capacitor `[MPa]`"
    p_leaf::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: define the struct to store the capacitor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the capacitor state and auxilary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ExtraXylemCapacitor{FT}
    "State variables of the capacitor"
    state::ExtraXylemCapacitorState{FT} = ExtraXylemCapacitorState{FT}()
    "Auxilary variables of the capacitor"
    auxil::ExtraXylemCapacitorAuxil{FT} = ExtraXylemCapacitorAuxil{FT}()
end
