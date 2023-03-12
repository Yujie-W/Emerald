module EmeraldEarth

using LazyArtifacts

using Dates: isleapyear
using Distributed: @everywhere, pmap
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using ProgressMeter: @showprogress

using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT
using GriddingMachine.Fetcher: ERA5SingleLevelsHourly, fetch_data!

using ..EmeraldCore.EarthGeometry: solar_zenith_angle
using ..EmeraldCore.Namespace: BetaFunction, BetaParameterG1, BetaParameterPsoil, MedlynSM, MonoMLTreeSPAC, Soil
using ..EmeraldCore.PhysicalChemistry: saturation_vapor_pressure
using ..EmeraldCore.SPAC: initialize!, soil_plant_air_continuum!, update!
using ..EmeraldIO.Netcdf: read_nc, save_nc!
using ..EmeraldIO.Text: read_csv
using ..EmeraldMath.Stats: nanmax, nanmean
using ..EmeraldUtility.Email: send_email!
using ..EmeraldUtility.Log: @tinfo
using ..EmeraldUtility.Threading: dynamic_workers!
using ..EmeraldUtility.Time: MDAYS, MDAYS_LEAP

# simulation settings
RESULT_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
SETUP_FOLDER  = "/home/wyujie/DATASERVER/model/CLIMA/LAND/setups";

# ERA5 settings
ERA5_FOLDER = "/home/wyujie/DATASERVER/reanalysis/ERA5/SingleLevels";


include("ear5.jl"           )
include("griddingmachine.jl")
include("grids.jl"          )
include("prescribe.jl"      )
include("simulation.jl"     )
include("threads.jl"        )


end # module
