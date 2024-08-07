# This file contains the structs used to store the energy state and auxiliary variables of the xylem (root and stem)
# Note here that the energy state needs to be initialized with the heat capacity of xylem and water, water content, and temperature...

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add XylemEnergyState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root energy state variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct XylemEnergyState{FT}
    "Total energy `[J]`"
    Σe::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2024-Feb-27: add XylemEnergySDAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root energy auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct XylemEnergySDAuxil{FT}
    "Combined heat capacity of root and water `[J K⁻¹]`"
    cp::FT = 0
    "Temperature `[K]`"
    t::FT = 298.15
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add XylemEnergyAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root energy auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct XylemEnergyAuxil{FT}
    "Partial derivative of the energy per time `[J s⁻¹]`"
    ∂e∂t::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add XylemEnergy
#     2024-Feb-27: add s_aux
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root energy variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct XylemEnergy{FT}
    "Xylem energy state"
    state::XylemEnergyState{FT} = XylemEnergyState{FT}()
    "State dependent auxiliary variables"
    s_aux::XylemEnergySDAuxil{FT} = XylemEnergySDAuxil{FT}()
    "Xylem energy auxil"
    auxil::XylemEnergyAuxil{FT} = XylemEnergyAuxil{FT}()
end;
