module EarthGeometry

using ..Constant: R_EQUATOR, R_POLAR, YEAR_D


include("earth.jl");
include("geometry.jl");
include("solar.jl");


end; # module
