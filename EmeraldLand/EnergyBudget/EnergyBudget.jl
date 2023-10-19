# This module is supposed to run the erngy budget of the SPAC (with and without plants)

module EnergyBudget

using ..EmeraldPhysics.Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, M_H₂O, ρ_H₂O

using ..Namespace: XylemHydraulics
using ..Namespace: JunctionCapacitor, Leaf, Root, Stem
using ..Namespace: AirLayer, SoilLayer
using ..Namespace: BulkSPAC

using ..PhysicalChemistry: latent_heat_vapor

using ..PlantHydraulics: flow_in, flow_out


include("cp.jl");
include("soil.jl");
include("junction.jl");
include("leaf.jl");
include("root.jl");
include("spac.jl");
include("stem.jl");


end; # module
