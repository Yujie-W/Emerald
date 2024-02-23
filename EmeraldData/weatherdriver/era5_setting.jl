# ERA5 settings
ERA5_FOLDER = "/home/wyujie/DATASERVER/reanalysis/ERA5/SingleLevels";
ERA5_LABELS = [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "2m_dewpoint_temperature",
            "2m_temperature",
            "mean_surface_direct_short_wave_radiation_flux",
            "mean_surface_direct_short_wave_radiation_flux_clear_sky",
            "mean_surface_downward_long_wave_radiation_flux",
            "mean_surface_downward_long_wave_radiation_flux_clear_sky",
            "mean_surface_downward_short_wave_radiation_flux",
            "mean_surface_downward_short_wave_radiation_flux_clear_sky",
            "mean_surface_downward_uv_radiation_flux",
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
            "msdrswrf", "msdrswrfcs", "msdwlwrf", "msdwlwrfcs","msdwswrf", "msdwswrfcs", "msdwuvrf",
            "skt", "stl1", "stl2", "stl3", "stl4", "sp", "tcc", "tp", "swvl1", "swvl2", "swvl3", "swvl4"];
ERA5_NETCDF = [
            "WIND_X", "WIND_Y", "T_DEW", "T_AIR",
            "RAD_DIR", "RAD_DIR_CS", "RAD_LW", "RAD_LW_CS", "RAD", "RAD_CS", "RAD_UV",
            "T_LEAF", "T_SOIL_1", "T_SOIL_2", "T_SOIL_3", "T_SOIL_4", "P_ATM", "CLOUD", "PRECIP", "SWC_1", "SWC_2", "SWC_3", "SWC_4"];
