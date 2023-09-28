# This file contains the capacitor outside the xylem, e.g., leaf capacitor and root-trunk junction capacitor

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: define the struct to store the state variables used in junction capacitor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state variables for junction capacitor

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct JunctionCapacitorState{FT}
    "Pressure volume curve of the capacitor"
    pv::AbstractPVCurve{FT} = ExponentialPVCurve{FT}()
    "Current volume of the capacitor `[mol m⁻²]`"
    v_storage::FT = 5000
    "Capacitor maximum volume per basal area or per leaf area `[mol m⁻²]`"
    v_max::FT = 5000
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: define the struct to store the auxilary variables used in junction capacitor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxilary variables for junction capacitor

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct JunctionCapacitorAuxil{FT}
    "Pressure of the capacitor `[MPa]`"
    pressure::FT = 0
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
Base.@kwdef mutable struct JunctionCapacitor{FT}
    "State variables of the capacitor"
    state::JunctionCapacitorState{FT} = JunctionCapacitorState{FT}()
    "Auxilary variables of the capacitor"
    auxil::JunctionCapacitorAuxil{FT} = JunctionCapacitorAuxil{FT}()
end
