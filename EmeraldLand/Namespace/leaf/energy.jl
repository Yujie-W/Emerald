# This file contains the structs used to store the energy state and auxiliary variables of the leaf
# Note here that the energy state needs to be initialized with the heat capacity of leaf and water, water content, and temperature...

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-25: add LeafEnergyState
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root energy state variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafEnergyState{FT}
    "Total energy `[J]`"
    Σe::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-25: add LeafEnergyAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root energy auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafEnergyAuxil{FT}
    #
    # these variables are to be initialized at sub steps
    #
    "Combined heat capacity of root and water `[J K⁻¹]`"
    cp::FT = 0
    "Temperature `[K]`"
    t::FT = T₂₅()
    "Partial derivative of the energy per time `[J s⁻¹]`"
    ∂e∂t::FT = 0
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-25: add LeafEnergy
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save root energy variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafEnergy{FT}
    "Leaf energy state"
    state::LeafEnergyState{FT} = LeafEnergyState{FT}()
    "Leaf energy auxil"
    auxil::LeafEnergyAuxil{FT} = LeafEnergyAuxil{FT}()
end;
