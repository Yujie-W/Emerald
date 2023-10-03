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
include("leaf/td.jl");

include("colimit.jl");
include("fluorescence.jl");
include("light_limited.jl");
include("model.jl");
include("product_limited.jl");
include("rubisco_limited.jl");
include("temperature.jl");


end # module
