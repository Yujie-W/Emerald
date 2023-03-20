module EmeraldFrontier

using Dates: isleapyear

using ..EmeraldData.ERA5: weather_driver_file
using ..EmeraldIO.Netcdf: read_nc
using ..EmeraldUtility.Time: month_days


# Netcdf settings for output
DF_VARIABLES  = ["F_H2O", "F_CO2", "F_GPP", "BETA", "SIF683", "SIF740", "SIF757", "SIF771", "RED", "BLUE", "NIR", "NDVI", "EVI", "NIRvI", "NIRvR", "PAR", "PPAR"];


include("driver.jl")


end # EmeraldFrontier
