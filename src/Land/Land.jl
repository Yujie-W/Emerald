module Land

using DataFrames: DataFrame, DataFrameRow
using GriddingMachine.Indexer: LandDatasetLabels, WeatherDriverLabels, grid_dict, grid_weather
using NetcdfIO: read_nc, save_nc!
using OrderedCollections: OrderedDict
using PkgUtility.EarthGeometry: solar_azimuth_angle, solar_zenith_angle
using PkgUtility.MathTools: nanmax, nanmean, nanmin, resample
using PkgUtility.PrettyDisplay: pretty_display!
using PkgUtility.UniversalConstants: M_H₂O, T₀, ρ_H₂O, energy_to_photon
using ProgressMeter: @showprogress
using Statistics: mean

using ..Namespace
using ..Namespace: BulkSPAC, SPACConfig
using ..PlantHydraulics: flow_out
using ..StomatalModels: read_β
# using ..SPAC: SAP_VOLUME
# using ..SPAC: CNPP, ET_SOIL, ET_VEGE, GPP, LATENT_HEAT, LEAF_PCI, K_PLANT, LONGWAVE_OUT, NET_LONGWAVE, NET_SHORTWAVE, OCS, SENSIBLE_HEAT, SHORTWAVE_OUT
# using ..SPAC: APAR, PAR, PPAR, ΦD_ΦN, ΦF_ΦP, ΣSIF, ΣSIF_CHL, ΣSIF_LEAF
# using ..SPAC: MODIS_BLUE, MODIS_EVI, MODIS_NDVI, MODIS_NIR, MODIS_NIRv, MODIS_NIRvR, MODIS_RED, OCO2_SIF759, OCO2_SIF770, TROPOMI_SIF683, TROPOMI_SIF740
using ..SPAC: dull_aux!, initialize_spac!, prescribe_air!, prescribe_soil!, prescribe_traits!, push_t_history!, soil_plant_air_continuum!, t_aux!


include("type/param-func-mapper.jl");

include("quantity/beta.jl");
include("quantity/co2.jl");
include("quantity/et.jl");
include("quantity/gpp.jl");
include("quantity/ocs.jl");
include("quantity/par.jl");

include("setting/saving.jl");
include("setting/setting.jl");

include("simulation/1-config.jl");
include("simulation/2-spac.jl");
include("simulation/3-driver.jl");
include("simulation/4-result.jl");
include("simulation/5-prescribe.jl");
include("simulation/6-simulation.jl");
include("simulation/7-save.jl");


end; # Land
