ERA5_FOLDER = "/home/wyujie/DATASERVER/reanalysis/ERA5/SingleLevels";


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-11: add the struct for ERA5 weather driver
#     2023-Mar-11: add field W_TOT for precipitation
#     2023-Mar-29: add field L_RAD for longwave radiation
#     2024-Mar-07: add thrid field to each tuple to store the variable name to use in a gridded dataframe file (for EmeraldFrontier and ClimaLandPRO)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for ERA5 Single Levels weather driver to use in EmeraldEarth

$(TYPEDFIELDS)

"""
Base.@kwdef struct ERA5SingleLevelsDriver
    "Total cloud cover"
    CLOUD::Tuple{String,String,String} = ("tcc", "total_cloud_cover", "CLOUD")
    "Clear sky longwave radiation"
    CS_LW::Tuple{String,String,String} = ("msdwlwrfcs", "mean_surface_downward_long_wave_radiation_flux_clear_sky", "RAD_LW_CS")
    "Clear sky direct shortwave radiation"
    CSDIR::Tuple{String,String,String} = ("msdrswrfcs", "mean_surface_direct_short_wave_radiation_flux_clear_sky", "RAD_DIR_CS")
    "Clear sky total shortwave radiation"
    CSRAD::Tuple{String,String,String} = ("msdwswrfcs", "mean_surface_downward_short_wave_radiation_flux_clear_sky", "RAD_CS")
    "Downward longwave radiation"
    L_RAD::Tuple{String,String,String} = ("msdwlwrf", "mean_surface_downward_long_wave_radiation_flux", "RAD_LW")
    "Atmospheric pressure"
    P_ATM::Tuple{String,String,String} = ("sp", "surface_pressure", "P_ATM")
    "Direct shortwave radiation"
    S_ALL::Tuple{String,String,String} = ("msdwswrf", "mean_surface_downward_short_wave_radiation_flux", "RAD")
    "Direct radiation"
    S_DIR::Tuple{String,String,String} = ("msdrswrf", "mean_surface_direct_short_wave_radiation_flux", "RAD_DIR")
    "Soil water content"
    SWC_1::Tuple{String,String,String} = ("swvl1", "volumetric_soil_water_layer_1", "SWC_1")
    "Soil water content"
    SWC_2::Tuple{String,String,String} = ("swvl2", "volumetric_soil_water_layer_2", "SWC_2")
    "Soil water content"
    SWC_3::Tuple{String,String,String} = ("swvl3", "volumetric_soil_water_layer_3", "SWC_3")
    "Soil water content"
    SWC_4::Tuple{String,String,String} = ("swvl4", "volumetric_soil_water_layer_4", "SWC_4")
    "Air temperature"
    T_AIR::Tuple{String,String,String} = ("t2m", "2m_temperature", "T_AIR")
    "Dew temperature"
    T_DEW::Tuple{String,String,String} = ("d2m", "2m_dewpoint_temperature", "T_DEW")
    "Soil temperature"
    T_S_1::Tuple{String,String,String} = ("stl1", "soil_temperature_level_1", "T_SOIL_1")
    "Soil temperature"
    T_S_2::Tuple{String,String,String} = ("stl2", "soil_temperature_level_2", "T_SOIL_2")
    "Soil temperature"
    T_S_3::Tuple{String,String,String} = ("stl3", "soil_temperature_level_3", "T_SOIL_3")
    "Soil temperature"
    T_S_4::Tuple{String,String,String} = ("stl4", "soil_temperature_level_4", "T_SOIL_4")
    "Skin temperature"
    T_SKN::Tuple{String,String,String} = ("skt", "skin_temperature", "T_LEAF")
    "Total UV radiation"
    UVRAD::Tuple{String,String,String} = ("msdwuvrf", "mean_surface_downward_uv_radiation_flux", "RAD_UV")
    "Total precipitation in m"
    W_TOT::Tuple{String,String,String} = ("tp", "total_precipitation", "PRECIP")
    "Wind speed"
    WINDU::Tuple{String,String,String} = ("u10", "10m_u_component_of_wind", "WIND_X")
    "Wind speed"
    WINDV::Tuple{String,String,String} = ("v10", "10m_v_component_of_wind", "WIND_Y")
end;


#######################################################################################################################################################################################################
#
# Changes to the functions
# General
#     2024-Mar-07: add function grid_file_path
#     2024-Mar-07: add function original_file_path
#     2024-Mar-07: add function reprocessed_file_path
#
#######################################################################################################################################################################################################
"""

    grid_file_path(gm_dict::Dict{String,Any})
Return the path of the weather driver file, given
- `gm_dict` GriddingMachine data dictionary

"""
function grid_file_path(gm_dict::Dict{String,Any})
    lat_ind = gm_dict["LAT_INDEX"];
    lon_ind = gm_dict["LON_INDEX"];
    nx      = gm_dict["RESO_SPACE"]
    year    = gm_dict["YEAR"];
    nc_name = "weather_driver_wd1_$(year)_$(lat_ind)_$(lon_ind)_$(nx)X.nc";

    return "$(DRIVER_FOLDER)/$(year)/$(nc_name)"
end;


"""

    original_file_path(gm_dict::Dict{String,Any}, varlabel::String)
    original_file_path(varlabel::String, year::Int)

Return the path of the original file, given
- `gm_dict` GriddingMachine data dictionary
- `varlabel` Variable label
- `year` Year

"""
function original_file_path end;

original_file_path(gm_dict::Dict{String,Any}, varlabel::String) = original_file_path(varlabel, gm_dict["YEAR"]);

original_file_path(varlabel::String, year::Int) = "$(ERA5_FOLDER)/original/$(varlabel)_SL_$(year).nc";


"""

    reprocessed_file_path(gm_dict::Dict{String,Any}, varlabel::String)
    reprocessed_file_path(varlabel::String, year::Int, nx::Int)

Return the path of the reprocessed file, given
- `gm_dict` GriddingMachine data dictionary
- `varlabel` Variable label
- `year` Year
- `nx` Number of grids in the 1 degree lat/lon

"""
function reprocessed_file_path end;

reprocessed_file_path(gm_dict::Dict{String,Any}, varlabel::String) = reprocessed_file_path(varlabel, gm_dict["YEAR"], gm_dict["RESO_SPACE"]);

reprocessed_file_path(varlabel::String, year::Int, nx::Int) = "$(ERA5_FOLDER)/reprocessed/$(varlabel)_SL_$(year)_$(nx)X.nc";
