using Test
import Emerald.EmeraldData.GlobalDatasets as GD
import Emerald.EmeraldFrontier as EF


@testset verbose = true "EmeraldFrontier" begin
    # Firstly, users need to prepare the files required first (download the files using provided functions in EmeraldData).
    # Secondly, users need to specify
    #     - Weather driver tag ("wd1" for ERA5 Single Levels data)
    #     - Dict of SPAC parameters (here from GriddingMachine)
    # Path to the weather driver data will be automatically retrieved with the tag and dict (which contains lat and lon information).
    gm_tag = "gm4";
    wd_tag = "wd1";
    gm_dict = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2019), 35.1, 115.2);
    df_simu = EF.simulation!(wd_tag, gm_dict; appending = false, selection = 1:24);
    @test true;
end;
