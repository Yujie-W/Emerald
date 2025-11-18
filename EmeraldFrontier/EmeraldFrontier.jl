module EmeraldFrontier

using DataFrames: DataFrame, DataFrameRow
using ProgressMeter: @showprogress
using Statistics: mean

using EmeraldUtilities.EarthGeometry: solar_azimuth_angle, solar_zenith_angle
using EmeraldUtilities.UniversalConstants: M_H₂O, T₀, ρ_H₂O
using NetcdfIO: read_nc, save_nc!

using ..EmeraldData.GlobalDatasets: grid_spac
using ..EmeraldData.WeatherDrivers: grid_weather_driver
using ..EmeraldLand.Namespace: BulkSPAC, SPACConfiguration
using ..EmeraldLand.SPAC: SAP_VOLUME
using ..EmeraldLand.SPAC: BETA, CNPP, ET_SOIL, ET_VEGE, GPP, LATENT_HEAT, K_PLANT, LONGWAVE_OUT, NET_LONGWAVE, NET_SHORTWAVE, OCS, SENSIBLE_HEAT, SHORTWAVE_OUT
using ..EmeraldLand.SPAC: APAR, PAR, PPAR, ΦD_ΦN, ΦF_ΦP, ΣSIF, ΣSIF_CHL, ΣSIF_LEAF
using ..EmeraldLand.SPAC: MODIS_BLUE, MODIS_EVI, MODIS_NDVI, MODIS_NIR, MODIS_NIRv, MODIS_NIRvR, MODIS_RED, OCO2_SIF759, OCO2_SIF770, TROPOMI_SIF683, TROPOMI_SIF740
using ..EmeraldLand.SPAC: dull_aux!, initialize_spac!, prescribe_air!, prescribe_soil!, prescribe_traits!, push_t_history!, soil_plant_air_continuum!, t_aux!
using ..EmeraldMath.Stats: nanmax, nanmean, nanmin


include("config.jl");
include("prepare_df.jl");
include("prescribe.jl");
include("simulation.jl");
include("save_fields.jl");


end; # EmeraldFrontier
