#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-11: add the struct for ERA5 weather driver
#     2023-Mar-11: add field W_TOT for precipitation
#     2023-Mar-29: add field L_RAD for longwave radiation
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for ERA5 Single Levels weather driver to use in EmeraldEarth

$(TYPEDFIELDS)

"""
Base.@kwdef struct ERA5SingleLevelsDriver
    "Downward longwave radiation"
    L_RAD::Tuple{String,String} = ("msdwlwrf", "mean_surface_downward_long_wave_radiation_flux")
    "Atmospheric pressure"
    P_ATM::Tuple{String,String} = ("sp", "surface_pressure")
    "Direct shortwave radiation"
    S_ALL::Tuple{String,String} = ("msdwswrf", "mean_surface_downward_short_wave_radiation_flux")
    "Direct radiation"
    S_DIR::Tuple{String,String} = ("msdrswrf", "mean_surface_direct_short_wave_radiation_flux")
    "Soil water content"
    SWC_1::Tuple{String,String} = ("swvl1", "volumetric_soil_water_layer_1")
    "Soil water content"
    SWC_2::Tuple{String,String} = ("swvl2", "volumetric_soil_water_layer_2")
    "Soil water content"
    SWC_3::Tuple{String,String} = ("swvl3", "volumetric_soil_water_layer_3")
    "Soil water content"
    SWC_4::Tuple{String,String} = ("swvl4", "volumetric_soil_water_layer_4")
    "Air temperature"
    T_AIR::Tuple{String,String} = ("t2m", "2m_temperature")
    "Dew temperature"
    T_DEW::Tuple{String,String} = ("d2m", "2m_dewpoint_temperature")
    "Soil temperature"
    T_S_1::Tuple{String,String} = ("stl1", "soil_temperature_level_1")
    "Soil temperature"
    T_S_2::Tuple{String,String} = ("stl2", "soil_temperature_level_2")
    "Soil temperature"
    T_S_3::Tuple{String,String} = ("stl3", "soil_temperature_level_3")
    "Soil temperature"
    T_S_4::Tuple{String,String} = ("stl4", "soil_temperature_level_4")
    "Skin temperature"
    T_SKN::Tuple{String,String} = ("skt", "skin_temperature")
    "Total UV radiation"
    UVRAD::Tuple{String,String} = ("msdwuvrf", "mean_surface_downward_uv_radiation_flux")
    "Total precipitation in m"
    W_TOT::Tuple{String,String} = ("tp", "total_precipitation")
    "Wind speed"
    WINDU::Tuple{String,String} = ("u10", "10m_u_component_of_wind")
    "Wind speed"
    WINDV::Tuple{String,String} = ("v10", "10m_v_component_of_wind")
end;
