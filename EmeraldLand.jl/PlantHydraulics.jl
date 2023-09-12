module PlantHydraulics

using Statistics: mean

using ..EmeraldMath.Solver: NewtonBisectionMethod, SolutionTolerance, find_zero
using ..EmeraldMath.Stats: nanmax, nanmean, nanmin

using ..Constant: CP_D_MOL, CP_L_MOL, CP_V_MOL, GAS_R, M_H₂O, T₂₅, ρg_MPa
using ..PhysicalChemistry: latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using ..Namespace: AbstractSoilVC, AbstractXylemVC, ComplexVC, LinearPVCurve, LogisticVC, PowerVC, SegmentedPVCurve, WeibullVC
using ..Namespace: BetaFunction, BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ
using ..Namespace: AbstractStomataModel, AndereggSM, BallBerrySM, EllerSM, GentineSM, LeuningSM, MedlynSM, SperrySM, WangSM, Wang2SM
using ..Namespace: Leaf, LeafHydraulics, Leaves1D, Leaves2D, NonSteadyStateFlow, Root, RootHydraulics, Soil, SoilLayer, SteadyStateFlow, Stem, StemHydraulics
using ..Namespace: MonoElementSPAC, MultiLayerSPAC, SPACConfiguration
using ..SoilHydraulics: soil_θ, soil_ψ_25

import ..SoilHydraulics: relative_hydraulic_conductance


# include functions related to beta factor
include("hydraulics/beta/model.jl");
include("hydraulics/beta/read.jl");
include("hydraulics/beta/set.jl");

# include functions related to drought effects
include("hydraulics/drought/disconnection.jl");
include("hydraulics/drought/legacy.jl");

# include function related to stomtal optimality
include("hydraulics/optimality/critical_flow.jl");
include("hydraulics/optimality/derivative.jl");

# include functions related to xylem vulnerability curve
include("hydraulics/vc/conductance.jl");
include("hydraulics/vc/pressure.jl");
include("hydraulics/vc/pressure_volume.jl");

# include functions related to flow and pressure profiles
include("hydraulics/flow_profile/flow_out.jl");
include("hydraulics/flow_profile/leaf_flow_out.jl");
include("hydraulics/flow_profile/leaf_flow_profile.jl");
include("hydraulics/flow_profile/read.jl");
include("hydraulics/flow_profile/root_pk.jl");
include("hydraulics/flow_profile/root_flow_out.jl");
include("hydraulics/flow_profile/root_flow_profile.jl");
include("hydraulics/flow_profile/spac_flow_profile.jl");
include("hydraulics/flow_profile/stem_flow_out.jl");
include("hydraulics/flow_profile/stem_flow_profile.jl");

include("hydraulics/budget.jl");
include("hydraulics/pressure_profile.jl");


end # module
