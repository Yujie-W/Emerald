module PlantHydraulics

using Statistics: mean

using ..EmeraldMath.Solver: NewtonBisectionMethod, SolutionTolerance, find_zero
using ..EmeraldMath.Stats: nanmax, nanmean, nanmin
using ..EmeraldMath.Math: lower_quadratic, upper_quadratic

using ..Constant: CP_D_MOL, CP_L_MOL, CP_V_MOL, GAS_R, M_H₂O, T₂₅, ρg_MPa
using ..PhysicalChemistry: latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using ..Namespace: AbstractSoilVC, AbstractXylemVC, ComplexVC, ExponentialPVCurve, LinearPVCurve, LogisticVC, PowerVC, SegmentedPVCurve, WeibullVC
using ..Namespace: XylemHydraulics, XylemHydraulicsAuxilNSS, XylemHydraulicsAuxilSS
using ..Namespace: BetaFunction, BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ
using ..Namespace: AbstractStomataModel, AndereggSM, BallBerrySM, EllerSM, GentineSM, LeuningSM, MedlynSM, SperrySM, WangSM, Wang2SM
using ..Namespace: Leaf, LeafHydraulics, Leaves1D, Leaves2D, NonSteadyStateFlow, Root, RootHydraulics, Soil, SoilLayer, SteadyStateFlow, Stem, StemHydraulics
using ..Namespace: MonoElementSPAC, MultiLayerSPAC, SPACConfiguration
using ..SoilHydraulics: soil_θ, soil_ψ_25


include("xylem/flow_profile.jl");
include("xylem/pv.jl");
include("xylem/vc.jl");


# include functions related to beta factor
include("beta/model.jl");
include("beta/read.jl");
include("beta/set.jl");

# include functions related to drought effects
include("drought/disconnection.jl");
include("drought/legacy.jl");

# include function related to stomtal optimality
include("optimality/critical_flow.jl");
include("optimality/derivative.jl");

# include functions related to xylem vulnerability curve
include("vc/pressure_volume.jl");

# include functions related to flow and pressure profiles
include("flow_profile/flow_out.jl");
include("flow_profile/leaf_flow_out.jl");
include("flow_profile/leaf_flow_profile.jl");
include("flow_profile/read.jl");
include("flow_profile/root_pk.jl");
include("flow_profile/root_flow_out.jl");
include("flow_profile/root_flow_profile.jl");
include("flow_profile/spac_flow_profile.jl");
include("flow_profile/stem_flow_out.jl");
include("flow_profile/stem_flow_profile.jl");

include("budget.jl");
include("pressure_profile.jl");


end # module
