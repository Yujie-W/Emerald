# This file contains functions to calculate the energy budgets of the xylem (root and stem)

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-29: add function heat_capacitance
#     2023-Oct-06: add method for SoilLayer
#     2023-Oct-07: account for the runoff water in the heat capacitance of the top soil
#     2023-Oct-09: add method for AirLayer
#
#######################################################################################################################################################################################################
"""

    heat_capacitance(organ)

Return the heat capacitance of the organ (xylem, root, stem, and leaf)

"""
function heat_capacitance end;

# soil
heat_capacitance(soiltr::SoilLayerTrait{FT}, soilta::SoilLayerTDAuxil{FT}, soilst::SoilLayerState{FT}; runoff::FT = FT(0)) where {FT} = (
    cp_gas = (soilst.ns[3] * CP_V_MOL(FT) + (soilst.ns[1] + soilst.ns[2] + soilst.ns[4] + soilst.ns[5]) * CP_D_MOL(FT)) / soilta.δz;

    # runoff in mol m⁻² s⁻¹, convert it to kg and then

    return soiltr.ρ * soiltr.cp + soilst.θ * ρ_H₂O(FT) * CP_L(FT) + cp_gas + runoff * CP_L_MOL(FT) / soilta.δz
);

heat_capacitance(soil::SoilLayer{FT}; runoff::FT = FT(0)) where {FT} = heat_capacitance(soil.trait, soil.t_aux, soil.state, runoff = runoff);

# plant
heat_capacitance(xylem::XylemHydraulics{FT}) where {FT} = (
    # The heat capaciatance of the xylem is the sum of the heat capaciatance of the water stored in the xylem and the heat capaciatance of the xylem itself
    return sum(xylem.state.v_storage) * CP_L_MOL(FT) + xylem.trait.cp * xylem.trait.area * xylem.trait.l
);

heat_capacitance(root::Root{FT}) where {FT} = heat_capacitance(root.xylem);

heat_capacitance(junc::JunctionCapacitor{FT}) where {FT} = (
    # the heat capaciatance of the junction is that of all water stored in the junction
    return junc.state.v_storage * CP_L_MOL(FT)
);

heat_capacitance(root::Stem{FT}) where {FT} = heat_capacitance(root.xylem);

heat_capacitance(capst::ExtraXylemCapacitorState{FT}, xyltr::XylemHydraulicsTrait{FT}, biotr::LeafBioTrait{FT}) where {FT} = (
    # leaf heat capaciatance is the sum of the heat capaciatance of the water stored in the leaf and the heat capaciatance of the leaf itself
    # here convert lma from g cm⁻² to kg m⁻² with the factor 10
    return (capst.v_storage * CP_L_MOL(FT) + xyltr.cp * biotr.lma * 10) * xyltr.area
);

heat_capacitance(leaf::Leaf{FT}) where {FT} = heat_capacitance(leaf.capacitor.state, leaf.xylem.trait, leaf.bio.trait);

# air
heat_capacitance(airst::AirLayerState{FT}) where {FT} = (airst.ns[1] + airst.ns[2] + airst.ns[4] + airst.ns[5]) * CP_D_MOL(FT) + airst.ns[3] * CP_V_MOL(FT);

heat_capacitance(air::AirLayer{FT}) where {FT} = heat_capacitance(air.state);
