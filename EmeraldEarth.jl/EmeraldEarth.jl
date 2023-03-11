module EmeraldEarth

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using ProgressMeter: @showprogress

using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT
using GriddingMachine.Fetcher: ERA5SingleLevelsHourly, fetch_data!

using ..EmeraldCore.Namespace: BetaFunction, BetaParameterG1, BetaParameterPsoil, MedlynSM, MonoMLTreeSPAC, Soil
using ..EmeraldCore.SPAC: initialize!, soil_plant_air_continuum!, update!
using ..EmeraldIO.Netcdf: read_nc, save_nc!
using ..EmeraldIO.Text: read_csv
using ..EmeraldMath.Stats: nanmax, nanmean
using ..EmeraldUtility.Email: send_email!
using ..EmeraldUtility.Log: @tinfo

# simulation settings
RESULT_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
SETUP_FOLDER  = "/home/wyujie/DATASERVER/model/CLIMA/LAND/setups";


include("ear5.jl"           )
include("griddingmachine.jl")
include("grids.jl"          )
include("weather.jl"        )


end # module
