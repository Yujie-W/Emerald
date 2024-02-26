#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-26: add method for LeafEnergyState-dependent auxiliary variables
#     2024-Feb-26: add method for Leaf
#
#######################################################################################################################################################################################################
s_aux!(capst::ExtraXylemCapacitorState{FT}, xylst::XylemHydraulicsState{FT}, biotr::LeafBioTrait{FT}, enest::LeafEnergyState{FT}, enesa::LeafEnergySDAuxil{FT}) where {FT} = (
    enesa.cp = heat_capacitance(capst, xylst, biotr);
    enesa.t = enest.Σe / enesa.cp;

    return nothing
);

s_aux!(leaf::Leaf{FT}) where {FT} = (
    s_aux!(leaf.capacitor.state, leaf.xylem.state, leaf.bio.trait, leaf.energy.state, leaf.energy.s_aux);

    return nothing
);
