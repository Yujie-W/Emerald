module StomatalModels

using ..EmeraldPhysics.Constant: CP_D_MOL, K_STEFAN, M_H₂O

using ..Namespace: AbstractSoilVC, AbstractXylemVC
using ..Namespace: AbstractStomataModel, AndereggSM, BallBerrySM, BetaParameterG1, BetaParameterVcmax, EllerSM, GentineSM, LeuningSM, MedlynSM, SperrySM, WangSM, Wang2SM
using ..Namespace: AirLayer, C4VJPModel, GCO₂Mode, HyperLeafBio, Leaves2D
using ..Namespace: MultiLayerSPAC
using ..Photosynthesis: leaf_photosynthesis!, ∂R∂T
using ..PhysicalChemistry: latent_heat_vapor, relative_diffusive_coefficient, relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using ..PlantHydraulics: ∂E∂P
using ..SoilHydraulics: relative_hydraulic_conductance


include("conductance.jl");
include("empirical.jl");
include("limits.jl");
include("optimality.jl");


end # module
