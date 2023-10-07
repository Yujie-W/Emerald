module SoilHydraulics

using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak
using ..EmeraldPhysics.Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V, CP_V_MOL, GAS_R, GRAVITY, M_H₂O, Λ_THERMAL_H₂O, ρ_H₂O, ρg_MPa

using ..Namespace: MultiLayerSPAC, Root, SPACConfiguration, VanGenuchten, XylemHydraulicsAuxilNSS, XylemHydraulicsAuxilSS
using ..PhysicalChemistry: diffusive_coefficient, latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure

import ..Namespace: BrooksCorey


# vc
include("vc/fitting.jl");
include("vc/vulnerability.jl");

# soil liquid water mass flow
include("liquid/infiltration.jl");
include("liquid/root.jl");

# soil trace gas diffusion
include("gas/diffusion.jl");
include("gas/volume.jl");


include("budget.jl");
include("diffusion.jl");
include("infiltration.jl");
include("volume.jl");


end # module
