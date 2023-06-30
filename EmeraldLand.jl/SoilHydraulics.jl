module SoilHydraulics

using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak

using ..Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, GRAVITY, M_H₂O, Λ_THERMAL_H₂O, ρ_H₂O, ρg_MPa
using ..Namespace: MultiLayerSPAC, NonSteadyStateFlow, Root, SPACConfiguration, SteadyStateFlow, VanGenuchten
using ..PhysicalChemistry: diffusive_coefficient, latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure

import ..Namespace: BrooksCorey


include("soil/budget.jl");
include("soil/constructor.jl");
include("soil/diffusion.jl");
include("soil/infiltration.jl");
include("soil/sink.jl");
include("soil/volume.jl");
include("soil/vulnerability.jl");


end # module
