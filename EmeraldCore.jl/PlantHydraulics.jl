module PlantHydraulics

using Statistics: mean

using ..EmeraldMath.Solver: NewtonBisectionMethod, SolutionTolerance, find_zero

using ..Constant: CP_D_MOL, CP_L_MOL, GAS_R, M_H₂O, T₂₅, ρg_MPa
using ..PhysicalChemistry: latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using ..Namespace: AbstractSoilVC, AbstractXylemVC, ComplexVC, LinearPVCurve, LogisticVC, PowerVC, SegmentedPVCurve, WeibullVC
using ..Namespace: BetaFunction, BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ
using ..Namespace: AbstractStomataModel, AndereggSM, BallBerrySM, EllerSM, GentineSM, LeuningSM, MedlynSM, SperrySM, WangSM, Wang2SM
using ..Namespace: Leaf, LeafHydraulics, Leaves1D, Leaves2D, NonSteadyStateFlow, Root, RootHydraulics, Soil, SoilLayer, SteadyStateFlow, Stem, StemHydraulics
using ..Namespace: MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
using ..Soil: soil_θ, soil_ψ_25

import ..Soil: relative_hydraulic_conductance


# include functions
include("hydraulics/beta.jl"             )
include("hydraulics/budget.jl"           )
include("hydraulics/critical_pressure.jl")
include("hydraulics/derivative.jl"       )
include("hydraulics/flow_profile.jl"     )
include("hydraulics/legacy.jl"           )
include("hydraulics/pressure_profile.jl" )
include("hydraulics/pressure_volume.jl"  )
include("hydraulics/target_flow.jl"      )
include("hydraulics/vulnerability.jl"    )


end # module