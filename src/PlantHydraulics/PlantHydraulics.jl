module PlantHydraulics

using Statistics: mean

using PkgUtility.MathTools: NewtonBisectionMethod, SolutionTolerance, find_zero, upper_quadratic
using PkgUtility.PhysicalChemistry: relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using PkgUtility.UniversalConstants: GAS_R, ρg_MPa

using ..Namespace: ComplexVC, LogisticVC, PowerVC, WeibullVC
using ..Namespace: ExponentialPVCurve, LinearPVCurve, SegmentedPVCurve
using ..Namespace: ExtraXylemCapacitor, ExtraXylemCapacitorAuxil, ExtraXylemCapacitorState, ExtraXylemCapacitorTrait
using ..Namespace: XylemHydraulics, XylemHydraulicsAuxilNSS, XylemHydraulicsAuxilSS, XylemHydraulicsState, XylemHydraulicsTrait
using ..Namespace: CanopyLayer, JunctionCapacitor, Leaf, Root, Stem
using ..Namespace: SoilLayer
using ..Namespace: BulkSPAC, SPACCache, SPACConfig

using ..SoilHydraulics: relative_soil_k


# xylem
include("xylem/conductance.jl");
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

# junction
include("junction/water_budget.jl");

# stem (dependent on xylem)
include("stem/flow_profile.jl");
include("stem/pressure_profile.jl");
include("stem/water_budget.jl");

# leaf (dependent on xylem)
include("leaf/capacitor.jl");
include("leaf/flow_profile.jl");
include("leaf/pressure_profile.jl");
include("leaf/water_budget.jl");

# plant (dependent on root, junction, stem, leaf)
include("plant/flow_profile.jl");
include("plant/pressure_profile.jl");
include("plant/water_budget.jl");


# include functions related to drought effects
include("drought/disconnection.jl");
include("drought/legacy.jl");


# include functions related to growth and recovery
include("growth/allocation.jl");
include("growth/newxylem.jl");
include("growth/recovery.jl");


end; # module
