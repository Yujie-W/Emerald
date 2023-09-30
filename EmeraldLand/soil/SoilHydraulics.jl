module SoilHydraulics

using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak

using ..Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V, CP_V_MOL, GAS_R, GRAVITY, M_H₂O, Λ_THERMAL_H₂O, ρ_H₂O, ρg_MPa
using ..Namespace: MultiLayerSPAC, Root, SPACConfiguration, VanGenuchten, XylemHydraulicsAuxilNSS, XylemHydraulicsAuxilSS
using ..PhysicalChemistry: diffusive_coefficient, latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure

import ..Namespace: BrooksCorey


include("budget.jl");
include("constructor.jl");
include("diffusion.jl");
include("infiltration.jl");
include("sink.jl");
include("volume.jl");
include("vulnerability.jl");


end # module
