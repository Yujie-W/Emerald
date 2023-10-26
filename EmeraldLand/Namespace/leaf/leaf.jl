# This file contains the fields within a leaf

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: add Leaf struct
#     2023-Oct-03: add fields photosystem, flux
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains the fields within a leaf

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Leaf{FT}
    "Leaf biophysics struct"
    bio::LeafBio{FT}
    "Extraxylary capacitor struct"
    capacitor::ExtraXylemCapacitor{FT} = ExtraXylemCapacitor{FT}()
    "Leaf energy struct"
    energy::LeafEnergy{FT} = LeafEnergy{FT}()
    "Leaf flux struct"
    flux::LeafFlux{FT}
    "Photosynthesis system struct"
    photosystem::LeafPhotosystem{FT} = LeafPhotosystem{FT}()
    "Leaf xylem hydraulics struct"
    xylem::XylemHydraulics{FT}
end;


"""

    Leaf(config::SPACConfiguration{FT}) where {FT}

Return the leaf struct with initialized energy states, given
- `config` `SPACConfiguration` type struct

"""
Leaf(config::SPACConfiguration{FT}) where {FT} = Leaf{FT}(bio = LeafBio(config), flux = LeafFlux(config), xylem = XylemHydraulics(config));
