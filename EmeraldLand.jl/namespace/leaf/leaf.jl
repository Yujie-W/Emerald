# This file contains the fields within a leaf

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Sep-26: add Leaf struct
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
    capacitor::ExtraXylemCapacitor{FT}
    "Leaf energy struct"
    energy::LeafEnergy{FT}
    "Leaf xylem hydraulics struct"
    xylem::XylemHydraulics{FT}
end


"""

    Leaf(config::SPACConfiguration{FT}) where {FT}

Return the leaf struct with initialized energy states, given
- `config` `SPACConfiguration` type struct

"""
Leaf(config::SPACConfiguration{FT}) where {FT} = (
    l_bio = HyperLeafBio(config);
    l_capacitor = ExtraXylemCapacitor{FT}(state = ExtraXylemCapacitorState{FT}(pv = SegmentedPVCurve{FT}(), v_max = 5));
    l_energy = LeafEnergy{FT}();
    l_xylem = XylemHydraulics(config);

    # now update the energy state of the leaf before returning the leaf struct
    l_xylem.state.cp = 1780;
    l_capacitor.state.v_storage = l_capacitor.state.v_max * l_xylem.state.area;
    l_energy.auxil.cp = l_capacitor.state.v_storage * CP_L_MOL(FT) + l_bio.state.lma * l_xylem.state.area * l_xylem.state.cp;
    l_energy.state.energy = l_energy.auxil.cp * l_energy.auxil.t;

    return Leaf{FT}(
                bio       = l_bio,
                capacitor = l_capacitor,
                energy    = l_energy,
                xylem     = l_xylem
    )
);
