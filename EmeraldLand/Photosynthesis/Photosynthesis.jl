module Photosynthesis

using DataFrames: DataFrame

using ..EmeraldMath.Math: lower_quadratic, upper_quadratic
using ..EmeraldPhysics.Constant: F_O₂, GAS_R

using ..Namespace: BetaFunction, BetaParameterG1, BetaParameterVcmax
using ..Namespace: AbstractStomataModel, BallBerrySM, GentineSM, LeuningSM, MedlynSM
using ..Namespace: Arrhenius, ArrheniusPeak, Q10, Q10Peak
using ..Namespace: MinimumColimit, QuadraticColimit, SerialColimit, SquareColimit
using ..Namespace: GCO₂Mode, PCO₂Mode
using ..Namespace: KNFluoscenceModel, QLFluoscenceModel
using ..Namespace: C3CytoState, C3CytoTrait, C3VJPState, C3VJPTrait, C4CLMTrait, C4VJPState, C4VJPTrait, LeafPhotosystemAuxil
using ..Namespace: LeafPhotosystem
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
