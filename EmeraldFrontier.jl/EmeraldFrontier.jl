module EmeraldFrontier

using DataFrames: DataFrameRow
using Dates: isleapyear
using DocStringExtensions: TYPEDEF, TYPEDFIELDS

using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT

using ..EmeraldCore.EarthGeometry: solar_zenith_angle
using ..EmeraldCore.Namespace: MultiLayerSPAC, Soil
using ..EmeraldCore.SPAC: initialize!, soil_plant_air_continuum!, update!
using ..EmeraldData.ERA5: weather_driver_file
using ..EmeraldIO.Netcdf: read_nc
using ..EmeraldIO.Text: read_csv
using ..EmeraldMath.Stats: nanmax, nanmean
using ..EmeraldUtility.Time: month_days


# Netcdf settings for output
DF_VARIABLES  = ["F_H2O", "F_CO2", "F_GPP", "BETA", "SIF683", "SIF740", "SIF757", "SIF771", "RED", "BLUE", "NIR", "NDVI", "EVI", "NIRvI", "NIRvR", "PAR", "PPAR"];


include("driver.jl"         )
include("griddingmachine.jl")
include("simulation.jl"     )
include("spac.jl"           )


end # EmeraldFrontier
