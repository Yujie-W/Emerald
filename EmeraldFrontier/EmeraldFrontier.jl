module EmeraldFrontier

using DataFrames: DataFrame, DataFrameRow
using ProgressMeter: @showprogress
using Statistics: mean

using NetcdfIO: read_nc, save_nc!

using ..EmeraldData.GlobalDatasets: grid_spac
using ..EmeraldData.WeatherDrivers: grid_weather_driver
using ..EmeraldLand.Namespace: BulkSPAC, SPACConfiguration
using ..EmeraldLand.SPAC: BETA, CNPP, GPP, PAR, PPAR, T_VEG, ΦF_ΦP, ΣSIF, ΣSIF_CHL, ΣSIF_LEAF
using ..EmeraldLand.SPAC: MODIS_BLUE, MODIS_EVI, MODIS_NDVI, MODIS_NIR, MODIS_NIRv, MODIS_NIRvR, MODIS_RED, OCO2_SIF759, OCO2_SIF770, TROPOMI_SIF683, TROPOMI_SIF740
using ..EmeraldLand.SPAC: dull_aux!, initialize_spac!, prescribe_air!, prescribe_soil!, prescribe_traits!, push_t_history!, soil_plant_air_continuum!, t_aux!
using ..EmeraldMath.Stats: nanmax, nanmean, nanmin
using ..EmeraldPhysics.Constant: M_H₂O, ρ_H₂O
using ..EmeraldPhysics.EarthGeometry: solar_azimuth_angle, solar_zenith_angle


# Netcdf settings for output
SAVING_DICT = Dict{String, Any}(
    # Modeled soil water content and temperature
            "MOD_SWC"    => true,
            "MOD_T_SOIL" => true,
    # Modeled leaf temperature
            "MOD_T_LEAF" => false,
            "MOD_T_MMM"  => true,
    # Modeled CO2 and H2O fluxes
            "BETA"       => false,
            "CNPP"       => true,
            "GPP"        => true,
            "ET_VEG"     => true,
    # SIF (default is false)
            "SIF683"     => false,
            "SIF740"     => true,
            "SIF757"     => false,
            "SIF771"     => false,
            "ΣSIF"       => false,
            "ΣSIF_CHL"   => false,
            "ΣSIF_LEAF"  => false,
            "ΦFΦP"       => false,
    # VI (default is false)
            "NDVI"       => false,
            "EVI"        => false,
            "NIRvI"      => false,
            "NIRvR"      => false,
            "PAR"        => false,
            "PPAR"       => false,
);


include("config.jl");
include("simulation.jl");


end; # EmeraldFrontier
