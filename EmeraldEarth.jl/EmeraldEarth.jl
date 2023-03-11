module EmeraldEarth

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using JLD2: save
using ProgressMeter: @showprogress

using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT
using GriddingMachine.Fetcher: ERA5SingleLevelsHourly, fetch_data!

using ..EmeraldIO.Netcdf: read_nc, save_nc!
using ..EmeraldIO.Text: read_csv
using ..EmeraldMath.Stats: nanmax, nanmean
using ..EmeraldUtility.Email: send_email!

# simulation settings
RESULT_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
SETUP_FOLDER  = "/home/wyujie/DATASERVER/model/CLIMA/LAND/setups";

# CLM5 settings
CLM5_PFTG = [0, 2.35, 2.35, 2.35, 4.12, 4.12, 4.45, 4.45, 4.45, 4.7, 4.7, 4.7, 2.22, 5.25, 1.62, 5.79, 5.79] .* sqrt(1000);
CLM5_PFTS = ["not_vegetated",
             "needleleaf_evergreen_temperate",
             "needleleaf_evergreen_boreal",
             "needleleaf_deciduous_boreal",
             "broadleaf_evergreen_tropical",
             "broadleaf_evergreen_temperate",
             "broadleaf_deciduous_tropical",
             "broadleaf_deciduous_temperate",
             "broadleaf_deciduous_boreal",
             "evergreen_shrub",
             "deciduous_temperate_shrub",
             "deciduous_boreal_shrub",
             "c3_arctic_grass",
             "c3_non-arctic_grass",
             "c4_grass",
             "c3_crop",
             "c3_irrigated"];

# ERA5 settings
ERA5_FOLDER = "/home/wyujie/DATASERVER/reanalysis/ERA5/SingleLevels";
ERA5_LABELS = [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "2m_dewpoint_temperature",
            "2m_temperature",
            "mean_surface_downward_long_wave_radiation_flux",
            "mean_surface_downward_long_wave_radiation_flux_clear_sky",
            "mean_surface_direct_short_wave_radiation_flux",
            "mean_surface_direct_short_wave_radiation_flux_clear_sky",
            "mean_surface_downward_short_wave_radiation_flux",
            "mean_surface_downward_short_wave_radiation_flux_clear_sky",
            "skin_temperature",
            "soil_temperature_level_1",
            "soil_temperature_level_2",
            "soil_temperature_level_3",
            "soil_temperature_level_4",
            "surface_pressure",
            "total_cloud_cover",
            "total_precipitation",
            "volumetric_soil_water_layer_1",
            "volumetric_soil_water_layer_2",
            "volumetric_soil_water_layer_3",
            "volumetric_soil_water_layer_4"];
ERA5_LAYERS = [
            "u10", "v10", "d2m", "t2m",
            "msdwlwrf", "msdwlwrfcs", "msdrswrf", "msdrswrfcs", "msdwswrf", "msdwswrfcs",
            "skt", "stl1", "stl2", "stl3", "stl4", "sp", "tcc", "tp", "swvl1", "swvl2", "swvl3", "swvl4"];
ERA5_NETCDF = [
            "WIND_X", "WIND_Y", "T_DEW", "T_AIR",
            "RAD_LW", "RAD_LW_CS", "RAD_DIR", "RAD_DIR_CS", "RAD", "RAD_CS",
            "T_LEAF", "T_SOIL_1", "T_SOIL_2", "T_SOIL_3", "T_SOIL_4", "P_ATM", "CLOUD", "PRECIP", "SWC_1", "SWC_2", "SWC_3", "SWC_4"];


include("ear5.jl"           )
include("griddingmachine.jl")


end # module
