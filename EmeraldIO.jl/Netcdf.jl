module Netcdf

using DataFrames: DataFrame
using NCDatasets: Dataset, defDim, defVar


# constants
const ATTR_LAT   = Dict("description" => "Latitude", "unit" => "°");
const ATTR_LON   = Dict("description" => "Longitude", "unit" => "°");
const ATTR_CYC   = Dict("description" => "Cycle index", "unit" => "-");
const ATTR_ABOUT = Dict("about" => "This is a file generated using Netcdf module of EmeraldIO.jl",
                        "notes" => "EmeraldIO.jl uses NCDatasets.jl to read and write NC files");


include("netcdf/append.jl");
include("netcdf/create.jl");
include("netcdf/grow.jl");
include("netcdf/info.jl");
include("netcdf/read.jl");
include("netcdf/save.jl");


end # module
