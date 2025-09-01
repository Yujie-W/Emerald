#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-11: add the struct for ERA5 weather driver
#     2023-Mar-11: add field W_TOT for precipitation
#     2023-Mar-29: add field L_RAD for longwave radiation
#     2024-Mar-07: add third field to each tuple to store the variable name to use in a gridded dataframe file (for EmeraldFrontier and ClimaLandPRO)
#     2025-Sep-01: remove unncessary fields
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for ERA5 Single Levels weather driver to use in EmeraldEarth

$(TYPEDFIELDS)

"""
Base.@kwdef struct ERA5SingleLevelsDriver
    "Downward longwave radiation"
    L_RAD::Tuple{String,String,String} = ("msdwlwrf", "mean_surface_downward_long_wave_radiation_flux", "RAD_LW")
    "Atmospheric pressure"
    P_ATM::Tuple{String,String,String} = ("sp", "surface_pressure", "P_ATM")
    "Direct shortwave radiation"
    S_ALL::Tuple{String,String,String} = ("msdwswrf", "mean_surface_downward_short_wave_radiation_flux", "RAD")
    "Direct radiation"
    S_DIR::Tuple{String,String,String} = ("msdrswrf", "mean_surface_direct_short_wave_radiation_flux", "RAD_DIR")
    "Air temperature"
    T_AIR::Tuple{String,String,String} = ("t2m", "2m_temperature", "T_AIR")
    "Dew temperature"
    T_DEW::Tuple{String,String,String} = ("d2m", "2m_dewpoint_temperature", "T_DEW")
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

    return "$(LAND_DRIVER)/$(year)/$(nc_name)"
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

original_file_path(varlabel::String, year::Int; folder::String = ERA5_SL_HOURLY) = "$(folder)/original/$(varlabel)_SL_$(year).nc";


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

reprocessed_file_path(varlabel::String, year::Int, nx::Int; folder::String = ERA5_SL_HOURLY) = "$(folder)/reprocessed/$(varlabel)_SL_$(year)_$(nx)X.nc";
