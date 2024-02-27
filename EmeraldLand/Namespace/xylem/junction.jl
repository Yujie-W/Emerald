# This file contains the capacitor outside the xylem, e.g., leaf capacitor and root-trunk junction capacitor

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-26: add struct JunctionCapacitorTrait
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the trait variables for junction capacitor

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct JunctionCapacitorTrait{FT}
    "Pressure volume curve of the capacitor"
    pv::Union{ExponentialPVCurve{FT}, LinearPVCurve{FT}, SegmentedPVCurve{FT}} = ExponentialPVCurve{FT}()
    "Capacitor maximum volume per basal area or per leaf area `[mol]`"
    v_max::FT = 5000
end;


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
    "Total energy storage of the capacitor `[J]`"
    Σe::FT = 5000 * CP_L_MOL() * T₂₅()
    "Current volume of the capacitor `[mol]`"
    v_storage::FT = 5000
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-27: add struct JunctionCapacitorSDAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the state dependent auxiliary variables for junction capacitor

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct JunctionCapacitorSDAuxil{FT}
    "Heat capacitance of the capacitor `[J K⁻¹]`"
    cp::FT = 0
    "Pressure of the capacitor `[MPa]`"
    pressure::FT = 0
    "Temperature of the capacitor `[K]`"
    t::FT = T₂₅()
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: define the struct to store the auxiliary variables used in junction capacitor
#     2023-Sep-30: add fields ∂e∂t and ∂w∂t
#     2023-Oct-02: add field cp
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the auxiliary variables for junction capacitor

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct JunctionCapacitorAuxil{FT}
    "Partial derivative of the energy per time `[J s⁻¹]`"
    ∂e∂t::FT = 0
    "Partial derivative of the water per time `[mol s⁻¹]`"
    ∂w∂t::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: define the struct to store the capacitor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the capacitor state and auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct JunctionCapacitor{FT}
    "Trait variables of the capacitor"
    trait::JunctionCapacitorTrait{FT} = JunctionCapacitorTrait{FT}()
    "State variables of the capacitor"
    state::JunctionCapacitorState{FT} = JunctionCapacitorState{FT}()
    "State dependent auxiliary variables of the capacitor"
    s_aux::JunctionCapacitorSDAuxil{FT} = JunctionCapacitorSDAuxil{FT}()
    "Auxilary variables of the capacitor"
    auxil::JunctionCapacitorAuxil{FT} = JunctionCapacitorAuxil{FT}()
end;
