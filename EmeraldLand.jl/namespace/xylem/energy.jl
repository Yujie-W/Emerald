# This file contains the structs used to store the energy state and auxiliary variables of the root
# Note here that the energy state needs to be initialized with the heat capacity of root and water, water content, and temperature...

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
    energy::FT = 0
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
    "Combined heat capacity of root and water `[J K⁻¹]`"
    cp::FT = 0
    "Temperature `[K]`"
    t::FT = 298.15
    "Partial derivative of the energy per time `[J s⁻¹]`"
    ∂e∂t::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add XylemEnergy
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
    "Xylem energy auxil"
    auxil::XylemEnergyAuxil{FT} = XylemEnergyAuxil{FT}()
end;
