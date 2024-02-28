module EmeraldFrontier

using DataFrames: DataFrame, DataFrameRow
using ProgressMeter: @showprogress

using NetcdfIO: read_nc, save_nc!

using ..EmeraldData.GlobalDatasets: grid_spac
using ..EmeraldData.WeatherDrivers: weather_driver_file
using ..EmeraldLand.Namespace: BulkSPAC, SPACConfiguration
using ..EmeraldLand.SPAC: BETA, CNPP, GPP, PAR, PPAR, T_VEG, ΦDFNP, ΣSIF, ΣSIF_CHL, ΣSIF_LEAF
using ..EmeraldLand.SPAC: MODIS_BLUE, MODIS_EVI, MODIS_NDVI, MODIS_NIR, MODIS_NIRv, MODIS_NIRvR, MODIS_RED, OCO2_SIF759, OCO2_SIF770, TROPOMI_SIF683, TROPOMI_SIF740
using ..EmeraldLand.SPAC: dull_aux!, initialize_spac!, prescribe_air!, prescribe_soil!, prescribe_traits!, soil_plant_air_continuum!, t_aux!
using ..EmeraldMath.Data: interpolate_data
using ..EmeraldMath.Stats: nanmax, nanmean, nanmin
using ..EmeraldPhysics.Constant: M_H₂O, ρ_H₂O
using ..EmeraldPhysics.EarthGeometry: solar_zenith_angle


# Netcdf settings for output
DF_SIMULATIONS = ["MOD_SWC_1", "MOD_SWC_2", "MOD_SWC_3", "MOD_SWC_4", "MOD_T_L_MAX", "MOD_T_L_MEAN", "MOD_T_L_MIN", "MOD_T_S_1", "MOD_T_S_2", "MOD_T_S_3", "MOD_T_S_4"];
DF_VARIABLES   = ["F_H2O", "F_CO2", "F_GPP", "BETA", "SIF683", "SIF740", "SIF757", "SIF771", "RED", "BLUE", "NIR", "NDVI", "EVI", "NIRvI", "NIRvR", "PAR", "PPAR", "ΦD", "ΦF", "ΦN", "ΦP", "ΣSIF",
                  "ΣSIF_CHL", "ΣSIF_LEAF"];


include("config.jl");
include("driver.jl");
include("simulation.jl");


end; # EmeraldFrontier
