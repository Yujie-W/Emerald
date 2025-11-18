module WeatherDrivers

using DataFrames: DataFrame
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using ProgressMeter: @showprogress

using EmeraldUtilities.MathTools: nanmean, resample_data
using EmeraldUtilities.PhysicalChemistry: saturation_vapor_pressure
using EmeraldUtilities.PrettyDisplay: pretty_display!
using GriddingMachine.Fetcher: fetch_data!
using NetcdfIO: append_nc!, read_nc, save_nc!, varname_nc

using ..EmeraldIO.Folders: ERA5_SL_HOURLY, LAND_DRIVER


# ERA5 settings and functions
include("era5_type.jl");

include("era5_grid.jl");
include("era5_load.jl");
include("era5_regrid.jl");

# parser and utility functions
include("parser.jl");
include("snapshot.jl");
include("verification.jl");


end; # module
