module EmeraldFrontier

using DataFrames: DataFrameRow
using Dates: isleapyear
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using ProgressMeter: @showprogress

using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT

using ..EmeraldData.ERA5: weather_driver_file
using ..EmeraldIO.Netcdf: read_nc, save_nc!
using ..EmeraldIO.Text: read_csv
using ..EmeraldLand.CanopyOptics: MODIS_BLUE, MODIS_EVI, MODIS_NDVI, MODIS_NIR, MODIS_NIRv, MODIS_NIRvR, MODIS_RED, OCO2_SIF759, OCO2_SIF770, TROPOMI_SIF683, TROPOMI_SIF740
using ..EmeraldLand.Constant: M_H₂O, ρ_H₂O
using ..EmeraldLand.EarthGeometry: solar_zenith_angle
using ..EmeraldLand.Namespace: MultiLayerSPAC, SPACConfiguration, Soil
using ..EmeraldLand.SPAC: BETA, CNPP, GPP, PPAR, T_VEG, initialize!, soil_plant_air_continuum!, update!
using ..EmeraldMath.Stats: nanmax, nanmean, nanmin
using ..EmeraldUtility.Time: month_days


# Netcdf settings for output
DF_SIMULATIONS = ["MOD_SWC_1", "MOD_SWC_2", "MOD_SWC_3", "MOD_SWC_4", "MOD_T_L_MAX", "MOD_T_L_MEAN", "MOD_T_L_MIN", "MOD_T_S_1", "MOD_T_S_2", "MOD_T_S_3", "MOD_T_S_4"];
DF_VARIABLES   = ["F_H2O", "F_CO2", "F_GPP", "BETA", "SIF683", "SIF740", "SIF757", "SIF771", "RED", "BLUE", "NIR", "NDVI", "EVI", "NIRvI", "NIRvR", "PAR", "PPAR"];


include("config.jl");
include("driver.jl");
include("griddingmachine.jl");
include("simulation.jl");
include("spac.jl");


end # EmeraldFrontier
