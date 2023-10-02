# This module is supposed to run the erngy budget of the SPAC (with and without plants)

module EnergyBudget

using ..EmeraldPhysics.Constant: CP_L_MOL, CP_V_MOL, M_H₂O

using ..Namespace: JunctionCapacitor, Leaf, Leaves2D, Root, Stem, XylemHydraulics
using ..Namespace: MultiLayerSPAC
using ..PhysicalChemistry: latent_heat_vapor
using ..PlantHydraulics: flow_in, flow_out


include("cp.jl");
include("junction.jl");
include("leaf.jl");
include("root.jl");
include("spac.jl");
include("stem.jl");


end # module
