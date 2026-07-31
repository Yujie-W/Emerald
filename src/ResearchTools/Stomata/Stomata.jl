module Stomata

using PkgUtility.PhysicalChemistry: saturation_vapor_pressure

using ..Namespace: SPACCache, SPACConfig
using ..Namespace: AirLayer, Leaf
using ..EPhotosynthesis: leaf_photosynthesis!
using ..PlantHydraulics: capacitance_volume, leaf_pressure_profile!, leaf_water_budget!, set_flow_profile!
using ..StomatalModels: stomatal_conductance!, ∂g∂t!
using ..SPAC: substep_aux!

using ..LeafLevelSetup: leaf_level_spac_cache


include("optimality.jl");
include("steady-state-solver.jl");


end # module
