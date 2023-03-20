module Photosynthesis

using ..EmeraldMath.Math: lower_quadratic, upper_quadratic

using ..Constant: F_O₂, GAS_R
using ..Namespace: CytochromeReactionCenter, VJPReactionCenter, VanDerTolFluorescenceModel
using ..Namespace: Arrhenius, ArrheniusPeak, C3VJPModel, C3CytochromeModel, C4VJPModel, GCO₂Mode, PCO₂Mode, Q10, Q10Peak
using ..Namespace: MinimumColimit, QuadraticColimit, SerialColimit, SquareColimit
using ..Namespace: AbstractStomataModel, BallBerrySM, BetaFunction, BetaParameterG1, BetaParameterVcmax, GentineSM, LeuningSM, MedlynSM
using ..Namespace: AirLayer, Leaf, Leaves1D, Leaves2D
using ..Namespace: MonoElementSPAC, MultiLayerSPAC


include("photosynthesis/colimit.jl"        )
include("photosynthesis/etr.jl"            )
include("photosynthesis/fluorescence.jl"   )
include("photosynthesis/light_limited.jl"  )
include("photosynthesis/model.jl"          )
include("photosynthesis/product_limited.jl")
include("photosynthesis/rubisco_limited.jl")
include("photosynthesis/temperature.jl"    )


end # module
