module PlantHydraulics

using Statistics: mean

using ..EmeraldMath.Solver: NewtonBisectionMethod, SolutionTolerance, find_zero
using ..EmeraldMath.Stats: nanmax, nanmean, nanmin
using ..EmeraldMath.Math: lower_quadratic, upper_quadratic

using ..Constant: CP_D_MOL, CP_L_MOL, CP_V_MOL, GAS_R, M_H₂O, T₂₅, ρg_MPa
using ..PhysicalChemistry: latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using ..Namespace: AbstractSoilVC, AbstractXylemVC, ComplexVC, ExponentialPVCurve, LinearPVCurve, LogisticVC, PowerVC, SegmentedPVCurve, WeibullVC
using ..Namespace: ExtraXylemCapacitor, ExtraXylemCapacitorAuxil, ExtraXylemCapacitorState, XylemHydraulics, XylemHydraulicsAuxilNSS, XylemHydraulicsAuxilSS, XylemHydraulicsState
using ..Namespace: Root, Root2, SoilLayer
using ..Namespace: BetaFunction, BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ
using ..Namespace: AbstractStomataModel, AndereggSM, BallBerrySM, EllerSM, GentineSM, LeuningSM, MedlynSM, SperrySM, WangSM, Wang2SM
using ..Namespace: Leaf, Leaves2D, Soil, Stem, Stem2
using ..Namespace: MultiLayerSPAC, SPACConfiguration
using ..SoilHydraulics: relative_hydraulic_conductance, soil_θ, soil_ψ_25


# xylem
include("xylem/critical_flow.jl");
include("xylem/flow_profile.jl");
include("xylem/pressure_profile.jl");
include("xylem/pv.jl");
include("xylem/vc.jl");
include("xylem/water_budget.jl");

# root (dependent on xylem)
include("root/flow_profile.jl");
include("root/pressure_profile.jl");
include("root/rhizosphere.jl");
include("root/water_budget.jl");

# stem (dependent on xylem)
include("stem/flow_profile.jl");
include("stem/pressure_profile.jl");
include("stem/water_budget.jl");

# leaf (dependent on xylem)
include("leaf/capacitor.jl");
include("leaf/flow_profile.jl");
include("leaf/pressure_profile.jl");
include("leaf/water_budget.jl");

# junction
include("junction/water_budget.jl");

# plant (dependent on root, junction, stem, leaf)
include("plant/flow_profile.jl");
include("plant/pressure_profile.jl");
include("plant/water_budget.jl");




# include functions related to beta factor
include("beta/model.jl");
include("beta/read.jl");
include("beta/set.jl");

# include functions related to drought effects
include("drought/disconnection.jl");
include("drought/legacy.jl");

# include function related to stomtal optimality
include("optimality/derivative.jl");

include("budget.jl");
include("pressure_profile.jl");


end # module
