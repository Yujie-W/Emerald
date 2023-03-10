module EmeraldCore

using ..EmeraldIO


include("Constant.jl" )
include("Namespace.jl")

include("EarthGeometry.jl"    )
include("Optics.jl"           )
include("PhysicalChemistry.jl")


end
