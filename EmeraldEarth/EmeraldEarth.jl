module EmeraldEarth

using LazyArtifacts

using DataFrames: DataFrame
using Dates: isleapyear
using Distributed: @everywhere, pmap
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using ProgressMeter: @showprogress

using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT
using NetcdfIO: append_nc!, create_nc!, grow_nc!, read_nc

using ..EmeraldData.WeatherDrivers: ERA5_FOLDER, ERA5SingleLevelsDriver
using ..EmeraldData.GlobalDatasets: LandDatasets, grid_dict, query_griddingmachine_data
using ..EmeraldIO.Text: read_csv
using ..EmeraldLand.Namespace: BetaFunction, BetaParameterG1, BetaParameterPsoil, MedlynSM, BulkSPAC, BulkSPACStates, SPACConfiguration, sync_state!
using ..EmeraldLand.PhysicalChemistry: saturation_vapor_pressure
using ..EmeraldLand.SPAC: GPP, PPAR, initialize_spac!, initialize_states!, prescribe_air!, prescribe_soil!, prescribe_traits!, soil_plant_air_continuum!
using ..EmeraldMath.Stats: nanmax, nanmean
using ..EmeraldPhysics.EarthGeometry: solar_zenith_angle
using ..EmeraldUtility.Log: @tinfo
using ..EmeraldUtility.Threading: dynamic_workers!
using ..EmeraldUtility.Time: MDAYS, MDAYS_LEAP


# simulation settings
RESULT_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
SETUP_FOLDER  = "/home/wyujie/DATASERVER/model/CLIMA/LAND/setups";
CACHE_CONFIG  = nothing;
CACHE_SPAC    = nothing;
CACHE_STATE   = nothing;


include("setup.jl");
include("threads.jl");
include("initialize.jl");


include("simulation.jl");


include("save.jl");


end; # module
