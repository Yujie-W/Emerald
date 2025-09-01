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
gm_dict["LAI"] = rand(12);
df = EDATA.WeatherDrivers.grid_weather_driver(wd_tag, gm_dict; appending = true);
df = EDATA.WeatherDrivers.grid_weather_driver(wd_tag, gm_dict, "/path/to/file/weather_driver_wd1_2001_774_1845_1X.nc");
"""
