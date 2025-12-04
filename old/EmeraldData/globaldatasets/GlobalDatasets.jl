module GlobalDatasets

using DataFrames: DataFrame
using Dates: isleapyear
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using Statistics: mean

using PkgUtility.EarthGeometry: solar_azimuth_angle, solar_zenith_angle
using PkgUtility.MathTools: gapfill_data!, nanmax, nanmean
using PkgUtility.PrettyDisplay: pretty_display!
using PkgUtility.TimeParser: MDAYS, MDAYS_LEAP
using GriddingMachine.Blender: regrid
using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT

using ..EmeraldIO.Text: read_csv
using ..EmeraldLand.Namespace: BetaFunction, BetaParameterG1, BetaParameterPsoil, MedlynSM, BulkSPAC, BulkSPACStates, SPACConfiguration, WangSM
using ..EmeraldLand.SPAC: initialize_spac!, prescribe_air!, prescribe_soil!, prescribe_traits!


CCS_1Y = read_csv("$(@__DIR__)/../../data/CO2-1Y.csv");
CCS_1M = read_csv("$(@__DIR__)/../../data/CO2-1M.csv");


include("co2.jl");
include("clm.jl");
include("land_datasets.jl");
include("query_data.jl");

include("extend_data.jl");
include("grid_dict.jl");
include("grid_spac.jl");


end # module
