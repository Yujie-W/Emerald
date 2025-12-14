import Emerald.EmeraldData.GlobalDatasets as GD
import Emerald.EmeraldFrontier as EF
import Emerald.EmeraldData as EDATA
gm_tag = "gm2";
wd_tag = "wd1";

lat = 49.47644;
lon = -126.25;


gmd = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2001), lat, lon);
gmd["MESSAGE_LEVEL"] = 1;
df = EDATA.WeatherDrivers.grid_weather_driver(wd_tag, gmd);
df_result = EF.simulation!(wd_tag, gmd);





gmd["LAI"] = rand(12);
df2 = EDATA.WeatherDrivers.grid_weather_driver(wd_tag, gmd);


gmd = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2001), 30.82, 104.1756; verification=false);
gmd["VCMAX25"] = 50.0;
gmd["MESSAGE_LEVEL"] = 1;
df = EDATA.WeatherDrivers.grid_weather_driver(wd_tag, gmd);
df_result = EF.simulation!(wd_tag, gmd);


gmd = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2001), 43.82, 87.61; verification=false);
gmd["VCMAX25"] = 50.0;
gmd["MESSAGE_LEVEL"] = 1;
df = EDATA.WeatherDrivers.grid_weather_driver(wd_tag, gmd);
df_result = EF.simulation!(wd_tag, gmd);
