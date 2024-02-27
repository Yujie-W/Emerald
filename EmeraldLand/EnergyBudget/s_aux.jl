#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-26: add method for LeafEnergyState-dependent auxiliary variables
#
#######################################################################################################################################################################################################
s_aux!(capst::ExtraXylemCapacitorState{FT}, xyltr::XylemHydraulicsTrait{FT}, biotr::LeafBioTrait{FT}, enest::LeafEnergyState{FT}, enesa::LeafEnergySDAuxil{FT}) where {FT} = (
    enesa.cp = heat_capacitance(capst, xyltr, biotr);
    enesa.t = enest.Σe / enesa.cp;

    return nothing
);
