module StomatalModels

using ..EmeraldPhysics.Constant: CP_D_MOL, K_STEFAN, M_H₂O

using ..Namespace: AbstractSoilVC
using ..Namespace: AbstractXylemVC
using ..Namespace: BetaFunction, BetaParameterG1, BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterVcmax, BetaParameterΘ
using ..Namespace: AbstractStomataModel, AndereggSM, BallBerrySM, EllerSM, GentineSM, LeuningSM, MedlynSM, SperrySM, WangSM, Wang2SM
using ..Namespace: CanopyLayer, Leaf, Root
using ..Namespace: AirLayer, SoilLayer
using ..Namespace: BulkSPAC, SPACCache

using ..Photosynthesis: photosynthesis_only!, ∂R∂T
using ..PhysicalChemistry: latent_heat_vapor, relative_diffusive_coefficient, relative_surface_tension, saturation_vapor_pressure
using ..PlantHydraulics: flow_in, flow_out, relative_xylem_k, xylem_end_pressure


# include functions related to beta factor
include("beta/model.jl");
include("beta/read.jl");
include("beta/set.jl");


# empirical models
include("empirical/ballberry.jl");
include("empirical/gentine.jl");
include("empirical/leuning.jl");
include("empirical/medlyn.jl");


# optimality models
include("optimality/dade.jl");
include("optimality/dedp.jl");

include("optimality/anderegg.jl");
include("optimality/eller.jl");
include("optimality/sperry.jl");
include("optimality/wang.jl");
include("optimality/wang2.jl");


# nighttime stomatal conductance
include("nighttime/drde.jl");
include("nighttime/dtde.jl");
include("nighttime/wang.jl");


# prognostic models
include("prognostic/daytime.jl");
include("prognostic/nighttime.jl");

include("prognostic/limits.jl");
include("prognostic/profile.jl");

include("prognostic/budget.jl");



end; # module
