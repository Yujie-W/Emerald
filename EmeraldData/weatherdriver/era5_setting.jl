# ERA5 settings
ERA5_FOLDER = "/home/wyujie/DATASERVER/reanalysis/ERA5/SingleLevels";
ERA5_DATASETS = [
            "10m_u_component_of_wind_SL"                                    "u10"           "WIND_X";
            "10m_v_component_of_wind_SL"                                    "v10"           "WIND_Y";
            "2m_dewpoint_temperature_SL"                                    "d2m"           "T_DEW";
            "2m_temperature_SL"                                             "t2m"           "T_AIR";
            "mean_surface_downward_long_wave_radiation_flux_SL"             "msdwlwrf"      "RAD_LW";
            "mean_surface_downward_long_wave_radiation_flux_clear_sky_SL"   "msdwlwrfcs"    "RAD_LW_CS";
            "mean_surface_direct_short_wave_radiation_flux_SL"              "msdrswrf"      "RAD_DIR";
            "mean_surface_direct_short_wave_radiation_flux_clear_sky_SL"    "msdrswrfcs"    "RAD_DIR_CS";
            "mean_surface_downward_short_wave_radiation_flux_SL"            "msdwswrf"      "RAD";
            "mean_surface_downward_short_wave_radiation_flux_clear_sky_SL"  "msdwswrfcs"    "RAD_CS";
            "mean_surface_downward_uv_radiation_flux_SL"                    "msdwuvrf"      "RAD_UV";
            "skin_temperature_SL"                                           "skt"           "T_LEAF";
            "soil_temperature_level_1_SL"                                   "stl1"          "T_SOIL_1";
            "soil_temperature_level_2_SL"                                   "stl2"          "T_SOIL_2";
            "soil_temperature_level_3_SL"                                   "stl3"          "T_SOIL_3";
            "soil_temperature_level_4_SL"                                   "stl4"          "T_SOIL_4";
            "surface_pressure_SL"                                           "sp"            "P_ATM";
            "total_cloud_cover_SL"                                          "tcc"           "CLOUD";
            "total_precipitation_SL"                                        "tp"            "PRECIP";
            "volumetric_soil_water_layer_1_SL"                              "swvl1"         "SWC_1";
            "volumetric_soil_water_layer_2_SL"                              "swvl2"         "SWC_2";
            "volumetric_soil_water_layer_3_SL"                              "swvl3"         "SWC_3";
            "volumetric_soil_water_layer_4_SL"                              "swvl4"         "SWC_4";]
ERA5_LABELS = ERA5_DATASETS[:, 1];
ERA5_LAYERS = ERA5_DATASETS[:, 2];
ERA5_NETCDF = ERA5_DATASETS[:, 3];
