module WeatherDrivers

using DataFrames: DataFrame
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using ProgressMeter: @showprogress

using GriddingMachine.Fetcher: fetch_data!
using NetcdfIO: append_nc!, read_nc, save_nc!, varname_nc

using ..EmeraldLand.PhysicalChemistry: saturation_vapor_pressure
using ..EmeraldMath.Stats: nanmean
using ..EmeraldUtility.Log: @tinfo
using ..EmeraldUtility.Email: send_email!


# CliMA Land settings
DRIVER_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/drivers";

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
