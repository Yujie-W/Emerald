module EmeraldCore

using ..EmeraldIO
using ..EmeraldMath


include("Constant.jl" )
include("Namespace.jl")

include("EarthGeometry.jl"    )
include("Optics.jl"           )
include("PhysicalChemistry.jl")

include("CanopyOptics.jl")
include("LeafOptics.jl"  )
include("Soil.jl"        )


end
