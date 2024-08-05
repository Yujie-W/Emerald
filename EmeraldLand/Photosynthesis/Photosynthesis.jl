module Photosynthesis

using DataFrames: DataFrame
using Statistics: mean

using ..EmeraldMath.Math: lower_quadratic, upper_quadratic
using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak
using ..EmeraldMath.Stats: nanmin, rmse
using ..EmeraldPhysics.Constant: F_O₂, GAS_R

using ..Namespace: BetaFunction, BetaParameterG1, BetaParameterVcmax
using ..Namespace: AbstractStomataModel, BallBerrySM, GentineSM, LeuningSM, MedlynSM
using ..Namespace: Arrhenius, ArrheniusPeak, Q10, Q10Peak, Q10PeakHT, Q10PeakLTHT
using ..Namespace: MinimumColimit, QuadraticColimit, SerialColimit, SquareColimit
using ..Namespace: AcMethodC3VcmaxPi, AcMethodC4Vcmax
using ..Namespace: AjMethodC3JmaxPi, AjMethodC3VqmaxPi, AjMethodC4JPSII
using ..Namespace: ApMethodC3Inf, ApMethodC3Vcmax, ApMethodC4VcmaxPi, ApMethodC4VpmaxPi
using ..Namespace: CytochromeFluorescenceModel, KNFluorescenceModel, QLFluorescenceModel, QLFluorescenceModelHan
using ..Namespace: GeneralC3Trait, C3State
using ..Namespace: GeneralC4Trait, C4State
using ..Namespace: CanopyLayerPhotosystem, CanopyLayerPhotosystemAuxil, LeafPhotosystem, LeafPhotosystemAuxil
using ..Namespace: CanopyLayer, Leaf
using ..Namespace: AirLayer
using ..Namespace: BulkSPAC, SPACCache, SPACConfiguration

# these are for the fitting of the A-Ci curve
using ..Namespace: ColimitCJCLMC3, ColimitCJCLMC4, ColimitIPCLM, ColimitJCLM, ηCTDJohnson, ηLTDJohnson


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
include("plant/layer.jl");
include("plant/leaf.jl");
include("plant/plant.jl");

# function to fit the traits
include("fitting/aci.jl");


end; # module
