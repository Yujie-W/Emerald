module Photosynthesis

using ..EmeraldMath.Math: lower_quadratic, upper_quadratic
using ..EmeraldPhysics.Constant: F_O₂, GAS_R

using ..Namespace: Arrhenius, ArrheniusPeak, Q10, Q10Peak
using ..Namespace: GCO₂Mode, PCO₂Mode
using ..Namespace: C3Cyto, C3VJP, C4VJP
using ..Namespace: MinimumColimit, QuadraticColimit, SerialColimit, SquareColimit
using ..Namespace: AbstractStomataModel, BallBerrySM, BetaFunction, BetaParameterG1, BetaParameterVcmax, GentineSM, LeuningSM, MedlynSM
using ..Namespace: AirLayer, Leaves2D
using ..Namespace: MultiLayerSPAC


# leaf level model
include("leaf/etr.jl");
include("leaf/light_limited.jl");
include("leaf/rubisco_limited.jl");
include("leaf/product_limited.jl");
include("leaf/td.jl");

include("colimit.jl");
include("fluorescence.jl");
include("model.jl");
include("temperature.jl");


end # module
