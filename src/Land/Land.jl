module Land

using DataFrames: DataFrame, DataFrameRow
using GriddingMachine.Indexer: LandDatasetLabels, WeatherDriverLabels, grid_dict, grid_weather
using NetcdfIO: read_nc, save_nc!
using OrderedCollections: OrderedDict
using PkgUtility.EarthGeometry: solar_azimuth_angle, solar_zenith_angle
using PkgUtility.MathTools: interpolate_data, nanmax, nanmean, nanmin, read_spectrum, resample
using PkgUtility.PrettyDisplay: pretty_display!
using PkgUtility.UniversalConstants: M_H₂O, K_STEFAN, T₀, ρ_H₂O, energy_to_photon
using ProgressMeter: @showprogress
using Statistics: mean

using ..Namespace
using ..Namespace: BulkSPAC, MultiLayerCanopy, ReferenceSpectra, ShortwaveRadiation, SPACConfig
using ..PlantHydraulics: flow_out
using ..StomatalModels: read_β
using ..SPAC: dull_aux!, initialize_spac!, prescribe_air!, prescribe_soil!, prescribe_traits!, push_t_history!, soil_plant_air_continuum!, t_aux!


include("type/param-func-mapper.jl");

include("quantity/flux.jl");
include("quantity/heat.jl");
include("quantity/PAR.jl");
include("quantity/SIF.jl");
include("quantity/stomata.jl");
include("quantity/VI.jl");
include("quantity/yield.jl");

include("setting/saving.jl");
include("setting/setting.jl");

include("simulation/1-config.jl");
include("simulation/2-spac.jl");
include("simulation/3-driver.jl");
include("simulation/4-result.jl");
include("simulation/5-prescribe.jl");
include("simulation/6-simulation.jl");


end; # Land
