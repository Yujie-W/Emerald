module EmeraldEarth

using LazyArtifacts

using Dates: isleapyear
using Distributed: @everywhere, pmap
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using ProgressMeter: @showprogress

using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT

using ..EmeraldCore.EarthGeometry: solar_zenith_angle
using ..EmeraldCore.Namespace: BetaFunction, BetaParameterG1, BetaParameterPsoil, MedlynSM, MultiLayerSPAC, MultiLayerSPACConfiguration, MultiLayerSPACState, Soil
using ..EmeraldCore.PhysicalChemistry: saturation_vapor_pressure
using ..EmeraldCore.SPAC: GPP, PPAR, initialize!, soil_plant_air_continuum!, spac_state!, update!
using ..EmeraldData.ERA5: ERA5_FOLDER, ERA5SingleLevelsDriver
using ..EmeraldIO.Netcdf: append_nc!, create_nc!, grow_nc!, read_nc
using ..EmeraldIO.Text: read_csv
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


include("griddingmachine.jl")
include("driver.jl"         )
include("cache.jl"          )
include("save.jl"           )
include("simulation.jl"     )
include("threads.jl"        )


end # module
