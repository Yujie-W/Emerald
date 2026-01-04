module ACi

using DataFrames: DataFrame
using PkgUtility.MathTools: ReduceStepMethodND, SolutionToleranceND, find_peak, nanmax, nanmin, rmse

using ..Namespace: SPACConfig
using ..Namespace: LeafPhotosystem
using ..Namespace: GeneralC3Trait, GeneralC4Trait
using ..Namespace: AcMethodC3VcmaxPi, AcMethodC4Vcmax, AjMethodC3JmaxPi, AjMethodC3VqmaxPi, AjMethodC4JPSII, ApMethodC3Vcmax, ApMethodC4VcmaxPi, ApMethodC4VpmaxPi
using ..Namespace: AirLayer
using ..Photosynthesis: photosynthesis!, photosystem_temperature_dependence!, temperature_correction


include("curve.jl");
include("rmse.jl");
include("fit.jl");

include("procedure/c3-vcmax-jmax-vcmax.jl");
include("procedure/c3-vcmax-vqmax-vcmax.jl");
include("procedure/c4-vcmax-jpsii-vcmaxpi.jl");
include("procedure/c4-vcmax-jpsii-vpmaxpi.jl");

include("pipeline.jl");


end # module
