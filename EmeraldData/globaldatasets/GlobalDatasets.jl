module GlobalDatasets

using DataFrames: DataFrame
using Dates: isleapyear
using DocStringExtensions: TYPEDEF, TYPEDFIELDS

using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT

using ..EmeraldIO.Text: read_csv
using ..EmeraldLand.Namespace: BetaFunction, BetaParameterG1, BetaParameterPsoil, MedlynSM, BulkSPAC, BulkSPACStates, SPACConfiguration
using ..EmeraldLand.SPAC: initialize_spac!, prescribe_air!, prescribe_soil!, prescribe_traits!
using ..EmeraldMath.Data: interpolate_data!
using ..EmeraldMath.Stats: nanmax, nanmean
using ..EmeraldPhysics.EarthGeometry: solar_azimuth_angle, solar_zenith_angle
using ..EmeraldUtility.Log: @tinfo
using ..EmeraldUtility.Time: MDAYS, MDAYS_LEAP


CCS = read_csv("$(@__DIR__)/../../data/CO2-1Y.csv");


include("clm.jl");
include("land_datasets.jl");
include("query_data.jl");

include("extend_data.jl");
include("grid_dict.jl");
include("grid_spac.jl");


end # module
