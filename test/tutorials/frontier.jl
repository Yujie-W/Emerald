using Test
import Emerald.EmeraldData.GlobalDatasets as GD
import Emerald.EmeraldFrontier as EF
import Emerald.EmeraldData as EDATA


@testset verbose = true "EmeraldFrontier" begin
    # Firstly, users need to prepare the files required first (download the files using provided functions in EmeraldData).
    # Secondly, users need to specify
    #     - Weather driver tag ("wd1" for ERA5 Single Levels data)
    #     - Dict of SPAC parameters (here from GriddingMachine)
    # Path to the weather driver data will be automatically retrieved with the tag and dict (which contains lat and lon information).
    gm_tag = "gm2";
    wd_tag = "wd1";
    gm_dict = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2019), 38.74, -92.20);
    df_simu = EF.simulation!(wd_tag, gm_dict; appending = false, selection = 1:24);
    @test true;
end;


"""
import Emerald.EmeraldData.GlobalDatasets as GD
import Emerald.EmeraldFrontier as EF
import Emerald.EmeraldData as EDATA
gm_tag = "gm2";
wd_tag = "wd1";
gm_dict = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2001), 38.74, -92.20);
df1 = EDATA.WeatherDrivers.grid_weather_driver(wd_tag, gm_dict);

gm_dict["LAI"] = rand(12);
df2 = EDATA.WeatherDrivers.grid_weather_driver(wd_tag, gm_dict);

df = EDATA.WeatherDrivers.grid_weather_driver(wd_tag, gm_dict, "/path/to/file/weather_driver_wd1_2001_774_1845_1X.nc");

lats = [,,,]
lons = [,,,]
for (lat, lon) in zip(lats, lons)
    gm_dict = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2001), lat, lon);
    df = EF.simulation!(wd_tag, gm_dict);
end;


gm_dict = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2001), 30.82, 104.1756; verification=false);
gm_dict["VCMAX25"] = 50.0;
gm_dict["MESSAGE_LEVEL"] = 1;
df = EDATA.WeatherDrivers.grid_weather_driver(wd_tag, gm_dict);
df_result = EF.simulation!(wd_tag, gm_dict);

"""
