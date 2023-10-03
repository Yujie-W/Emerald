# This file contains the fields within a leaf

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: add Leaf struct
#     2023-Oct-03: add field photosystem
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
    bio::HyperLeafBio{FT}
    "Extraxylary capacitor struct"
    capacitor::ExtraXylemCapacitor{FT} = ExtraXylemCapacitor{FT}()
    "Leaf energy struct"
    energy::LeafEnergy{FT} = LeafEnergy{FT}()
    "Photosynthesis system struct"
    photosystem::Union{C3Cyto{FT}, C3VJP{FT}, C4VJP{FT}} = C3VJP{FT}()
    "Leaf xylem hydraulics struct"
    xylem::XylemHydraulics{FT}
end


"""

    Leaf(config::SPACConfiguration{FT}) where {FT}

Return the leaf struct with initialized energy states, given
- `config` `SPACConfiguration` type struct

"""
Leaf(config::SPACConfiguration{FT}) where {FT} = (
    leaf = Leaf{FT}(bio = HyperLeafBio(config), xylem = XylemHydraulics(config));
    initialize_energy_storage!(leaf);

    return leaf
);

initialize_energy_storage!(leaf::Leaf{FT}) where {FT} = (
    leaf.xylem.state.cp = 1780;
    leaf.xylem.state.k_max = 0.04;
    leaf.capacitor.state.v_storage = leaf.capacitor.state.v_max * leaf.xylem.state.area;
    leaf.energy.auxil.cp = leaf.capacitor.state.v_storage * CP_L_MOL(FT) + leaf.bio.state.lma * 10 * leaf.xylem.state.area * leaf.xylem.state.cp;
    leaf.energy.state.Σe = leaf.energy.auxil.cp * leaf.energy.auxil.t;

    return nothing
);
