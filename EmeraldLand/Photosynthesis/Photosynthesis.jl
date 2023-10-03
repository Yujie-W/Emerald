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


# photosystem level model (order by the step in photosynthesis model)
include("photosystem/td.jl");
include("photosystem/etr.jl");
include("photosystem/rubisco_limited.jl");
include("photosystem/light_limited.jl");
include("photosystem/product_limited.jl");
include("photosystem/colimit.jl");
include("photosystem/fluorescence.jl");

include("model.jl");
include("temperature.jl");


end # module
