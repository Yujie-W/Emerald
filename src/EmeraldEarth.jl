module EmeraldEarth

using LazyArtifacts

using Dates: isleapyear
using Distributed: @everywhere, pmap
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using ProgressMeter: @showprogress

using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT
using NetcdfIO: append_nc!, create_nc!, grow_nc!, read_nc

using ..EmeraldData.ERA5: ERA5_FOLDER, ERA5SingleLevelsDriver
using ..EmeraldIO.Text: read_csv
using ..EmeraldLand.EarthGeometry: solar_zenith_angle
using ..EmeraldLand.Namespace: BetaFunction, BetaParameterG1, BetaParameterPsoil, MedlynSM, MultiLayerSPAC, SPACConfiguration, MultiLayerSPACState, Soil
using ..EmeraldLand.PhysicalChemistry: saturation_vapor_pressure
using ..EmeraldLand.SPAC: GPP, PPAR, initialize!, soil_plant_air_continuum!, spac_state!, update!
using ..EmeraldMath.Stats: nanmax, nanmean
using ..EmeraldUtility.Log: @tinfo
using ..EmeraldUtility.Threading: dynamic_workers!
using ..EmeraldUtility.Time: MDAYS, MDAYS_LEAP


# simulation settings
RESULT_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
SETUP_FOLDER  = "/home/wyujie/DATASERVER/model/CLIMA/LAND/setups";
CACHE_CONFIG  = nothing;
CACHE_SPAC    = nothing;
CACHE_STATE   = nothing;


include("../EmeraldEarth.jl/griddingmachine.jl");
include("../EmeraldEarth.jl/driver.jl");
include("../EmeraldEarth.jl/cache.jl");
include("../EmeraldEarth.jl/save.jl");
include("../EmeraldEarth.jl/simulation.jl");
include("../EmeraldEarth.jl/threads.jl");


end # module