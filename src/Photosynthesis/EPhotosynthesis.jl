module EPhotosynthesis

using DataFrames: DataFrame
using Statistics: mean

using PkgUtility.MathTools: ReduceStepMethodND, SolutionToleranceND, find_peak, lower_quadratic, nanmax, nanmin, rmse, upper_quadratic
using PkgUtility.UniversalConstants: F_O₂, GAS_R

using Photosynthesis: Arrhenius, ArrheniusPeak, Q10, Q10Peak
using Photosynthesis: LeafPhotosystem
using Photosynthesis: colimit_photosynthesis!, light_limited_rate!, photosystem_coefficients!, photosystem_electron_transport!, photosystem_temperature_dependence!, product_limited_rate!, rubisco_limited_rate!, temperature_corrected_value

using ..Namespace: BetaFunction, BetaParameterG1, BetaParameterVcmax
using ..Namespace: AbstractStomatalConductanceModel, BallBerrySM, GentineSM, LeuningSM, MedlynSM
using ..Namespace: Leaf
using ..Namespace: AirLayer
using ..Namespace: BulkSPAC, SPACCache, SPACConfig

# these are for the fitting of the A-Ci curve
# using ..Namespace: Arrhenius, ArrheniusPeak, ArrheniusPeak2, Q10, Q10Peak, Q10PeakHT, Q10PeakLTHT
# using ..Namespace: MinimumColimit, QuadraticColimit, SerialColimit, SquareColimit, UnionColimit
# using ..Namespace: AcMethodC3VcmaxPi, AcMethodC4Vcmax
# using ..Namespace: AjMethodC3JmaxPi, AjMethodC3VqmaxPi, AjMethodC4JPSII
# using ..Namespace: ApMethodC3Inf, ApMethodC3Vcmax, ApMethodC4VcmaxPi, ApMethodC4VpmaxPi
# using ..Namespace: CytochromeFluorescenceModel, KNFluorescenceModel, QLFluorescenceModel, QLFluorescenceModelHan
# using ..Namespace: GeneralC3Trait, C3State
# using ..Namespace: GeneralC4Trait, C4State
# using ..Namespace: ColimitCJCLMC3, ColimitCJCLMC4, ColimitIPCLM, ColimitJCLM, ηCTDJohnson, ηLTDJohnson
# using ..Namespace: LeafPhotosystem, LeafPhotosystemAuxil

# functions to use with stomatal models
include("stomata/derivative.jl");
include("stomata/photo_only.jl");

# functions to use with SPAC
include("plant/layer.jl");
include("plant/plant.jl");


end; # module
