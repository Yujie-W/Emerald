module GlobalDatasets

using DataFrames: DataFrame
using DocStringExtensions: TYPEDEF, TYPEDFIELDS

using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT

using ..EmeraldIO.Text: read_csv
using ..EmeraldMath.Data: interpolate_data!
using ..EmeraldMath.Stats: nanmax, nanmean
using ..EmeraldUtility.Log: @tinfo


CCS = read_csv("$(@__DIR__)/../../data/CO2-1Y.csv");


include("clm.jl");
include("land_datasets.jl");

include("extend_data.jl");
include("grid_dict.jl");


end # module
