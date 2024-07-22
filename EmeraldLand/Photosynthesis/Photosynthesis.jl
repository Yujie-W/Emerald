module Photosynthesis

using DataFrames: DataFrame

using ..EmeraldMath.Math: lower_quadratic, upper_quadratic
using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak
using ..EmeraldMath.Stats: nanmin, rmse
using ..EmeraldPhysics.Constant: F_O₂, GAS_R

using ..Namespace: BetaFunction, BetaParameterG1, BetaParameterVcmax
using ..Namespace: AbstractStomataModel, BallBerrySM, GentineSM, LeuningSM, MedlynSM
using ..Namespace: Arrhenius, ArrheniusPeak, Q10, Q10Peak
using ..Namespace: MinimumColimit, QuadraticColimit, SerialColimit, SquareColimit
using ..Namespace: GCO₂Mode, PCO₂Mode
using ..Namespace: KNFluoscenceModel, QLFluoscenceModel
using ..Namespace: C3CytoTrait, C3CytoState
using ..Namespace: C3CLMTrait, C3FvCBTrait, C3VJPTrait, C3VJPState
using ..Namespace: C4CLMTrait, C4VJPTrait, C4VJPState
using ..Namespace: LeafPhotosystem, LeafPhotosystemAuxil
using ..Namespace: Leaf
using ..Namespace: AirLayer
using ..Namespace: BulkSPAC


# photosystem level model (order by the step in photosynthesis model)
include("photosystem/td.jl");
include("photosystem/etr.jl");
include("photosystem/rubisco_limited.jl");
include("photosystem/light_limited.jl");
include("photosystem/product_limited.jl");
include("photosystem/colimit.jl");
include("photosystem/fluorescence.jl");
include("photosystem/photosystem.jl");

# functions to use with stomatal models
include("stomata/derivative.jl");
include("stomata/photo_only.jl");

# functions to use with SPAC
include("plant/leaf.jl");
include("plant/plant.jl");

# function to fit the traits
include("fitting/aci.jl");


end; # module
