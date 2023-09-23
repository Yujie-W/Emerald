# This file contains the structs used to store the energy state and auxiliary variables of the root
# Note here that the energy state needs to be initialized with the heat capacity of root and water, water content, and temperature...

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add RootEnergyState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root energy state variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct RootEnergyState{FT}
    "Total energy `[J]`"
    energy::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add RootEnergyAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root energy auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct RootEnergyAuxil{FT}
    "Combined heat capacity of root and water `[J K⁻¹]`"
    cp::FT = 0
    "Temperature `[K]`"
    t::FT = 298.15
    "Partial derivative of the energy per time `[J s⁻¹]`"
    ∂e∂t::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-23: add RootEnergy
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root energy variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct RootEnergy{FT}
    "Root energy state"
    state::RootEnergyState{FT} = RootEnergyState{FT}()
    "Root energy auxil"
    auxil::RootEnergyAuxil{FT} = RootEnergyAuxil{FT}()
end
