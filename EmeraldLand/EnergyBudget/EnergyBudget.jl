# This module is supposed to run the erngy budget of the SPAC (with and without plants)

module EnergyBudget

using Statistics: mean

using ..EmeraldPhysics.Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, GLUCOSE, M_H₂O, ρ_H₂O

using ..Namespace: XylemHydraulics
using ..Namespace: CanopyLayer, JunctionCapacitor, Leaf, Root, Stem
using ..Namespace: AirLayer, AirLayerState
using ..Namespace: SoilLayer, SoilLayerState, SoilLayerTrait, SoilLayerTDAuxil
using ..Namespace: ExtraXylemCapacitorState, LeafBioTrait, LeafEnergyState, LeafEnergySDAuxil, XylemHydraulicsState, XylemHydraulicsTrait
using ..Namespace: BulkSPAC, SPACConfiguration
using ..PhysicalChemistry: latent_heat_vapor
using ..PlantHydraulics: flow_in, flow_out


include("cp.jl");
include("junction.jl");
include("leaf.jl");
include("root.jl");
include("soil.jl");
include("spac.jl");
include("stem.jl");


end; # module
