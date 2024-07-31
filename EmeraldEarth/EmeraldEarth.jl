module EmeraldEarth

using LazyArtifacts

using DataFrames: DataFrame
using Dates: isleapyear
using Distributed: @everywhere, pmap
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using ProgressMeter: @showprogress

using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT
using NetcdfIO: append_nc!, create_nc!, grow_nc!, read_nc

using ..EmeraldData.WeatherDrivers: ERA5SingleLevelsDriver, dict_contains_nan
using ..EmeraldData.GlobalDatasets: LandDatasets, grid_dict, grid_spac, prescribe_gm_wd_data!, query_griddingmachine_data
using ..EmeraldIO.Text: read_csv
using ..EmeraldLand.Namespace: BetaFunction, BetaParameterG1, BetaParameterPsoil, MedlynSM, BulkSPAC, BulkSPACStates, SPACConfiguration
using ..EmeraldLand.Namespace: sync_state!
using ..EmeraldLand.PhysicalChemistry: saturation_vapor_pressure
using ..EmeraldLand.SPAC: GPP, PPAR, initialize_spac!, prescribe_air!, prescribe_soil!, prescribe_traits!, push_t_history!, soil_plant_air_continuum!
using ..EmeraldMath.Stats: nanmax, nanmean
using ..EmeraldUtility.Log: @tinfo
using ..EmeraldUtility.Threading: dynamic_workers!
using ..EmeraldUtility.Time: MDAYS, MDAYS_LEAP


# simulation settings
RESULT_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
SETUP_FOLDER  = "/home/wyujie/DATASERVER/model/CLIMA/LAND/setups";
CACHE_CONFIG  = nothing;

include("setup.jl");
include("threads.jl");
include("initialize.jl");
include("simulation.jl");




include("save.jl");


end; # module
