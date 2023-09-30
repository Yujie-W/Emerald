module Photosynthesis

using ..EmeraldMath.Math: lower_quadratic, upper_quadratic

using ..Constant: F_O₂, GAS_R
using ..Namespace: CytochromeReactionCenter, VJPReactionCenter, VanDerTolFluorescenceModel
using ..Namespace: Arrhenius, ArrheniusPeak, C3VJPModel, C3CytochromeModel, C4VJPModel, GCO₂Mode, PCO₂Mode, Q10, Q10Peak
using ..Namespace: MinimumColimit, QuadraticColimit, SerialColimit, SquareColimit
using ..Namespace: AbstractStomataModel, BallBerrySM, BetaFunction, BetaParameterG1, BetaParameterVcmax, GentineSM, LeuningSM, MedlynSM
using ..Namespace: AirLayer, Leaves2D
using ..Namespace: MultiLayerSPAC


include("colimit.jl");
include("etr.jl");
include("fluorescence.jl");
include("light_limited.jl");
include("model.jl");
include("product_limited.jl");
include("rubisco_limited.jl");
include("temperature.jl");


end # module
