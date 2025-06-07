module SoilHydraulics

using ..EmeraldPhysics.Constant: CP_D_MOL, CP_I_MOL, CP_L, GAS_R, M_H₂O, T₀, ρ_H₂O, ρg_MPa

using ..Namespace: BrooksCorey, VanGenuchten
using ..Namespace: XylemHydraulicsAuxilNSS, XylemHydraulicsAuxilSS
using ..Namespace: Root
using ..Namespace: SoilLayer
using ..Namespace: BulkSPAC, SPACConfiguration
using ..PhysicalChemistry: diffusive_coefficient, latent_heat_melt, latent_heat_vapor, relative_surface_tension, saturation_vapor_pressure


# vc
include("vc/fitting.jl");
include("vc/vulnerability.jl");

# soil trace gas diffusion
include("gas/diffusion.jl");
include("gas/volume.jl");

# soil liquid water mass flow
include("liquid/infiltration.jl");
include("liquid/root.jl");
include("liquid/runoff.jl");

# soil flow profiles
include("soil/profile.jl");
include("soil/budget.jl");


end; # module
