module Land

using DataFrames: DataFrame, DataFrameRow
using GriddingMachine.Indexer: LandDatasetLabels, WeatherDriverLabels, grid_dict, grid_weather
using NetcdfIO: read_nc, save_nc!
using PkgUtility.EarthGeometry: solar_azimuth_angle, solar_zenith_angle
using PkgUtility.MathTools: nanmax, nanmean, nanmin, resample
using PkgUtility.UniversalConstants: M_H₂O, T₀, ρ_H₂O
using ProgressMeter: @showprogress
using Statistics: mean

using ..Namespace: BulkSPAC, SPACConfig
using ..StomatalModels: WangSM
using ..SPAC: SAP_VOLUME
using ..SPAC: BETA, CNPP, ET_SOIL, ET_VEGE, GPP, LATENT_HEAT, LEAF_PCI, K_PLANT, LONGWAVE_OUT, NET_LONGWAVE, NET_SHORTWAVE, OCS, SENSIBLE_HEAT, SHORTWAVE_OUT
using ..SPAC: APAR, PAR, PPAR, ΦD_ΦN, ΦF_ΦP, ΣSIF, ΣSIF_CHL, ΣSIF_LEAF
using ..SPAC: MODIS_BLUE, MODIS_EVI, MODIS_NDVI, MODIS_NIR, MODIS_NIRv, MODIS_NIRvR, MODIS_RED, OCO2_SIF759, OCO2_SIF770, TROPOMI_SIF683, TROPOMI_SIF740
using ..SPAC: dull_aux!, initialize_spac!, prescribe_air!, prescribe_soil!, prescribe_traits!, push_t_history!, soil_plant_air_continuum!, t_aux!


include("initialize-config.jl");
include("initialize-driver.jl");
include("initialize-saving.jl");
include("initialize-spac.jl");

include("site-prescribe.jl");
include("site-simulation.jl");

include("saving-dict.jl");
include("saving-parser.jl");


end; # Land
