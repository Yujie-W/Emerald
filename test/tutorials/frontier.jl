@testset verbose = true "EmeraldFrontier" begin
    # Firstly, users need to prepare the files required first (download the files using provided functions in EmeraldData).
    # Secondly, users need to specify
    #     - Weather driver tag ("wd1" for ERA5 Single Levels data)
    #     - Dict of SPAC parameters (here from GriddingMachine)
    # Path to the weather driver data will be automatically retrieved with the tag and dict (which contains lat and lon information).
    wd_tag = "wd1";
    gm_dict = EmeraldFrontier.gm_dict(EmeraldFrontier.GriddingMachineLabels(year=2019), 35.1, 115.2);
    wd_data = EmeraldFrontier.weather_driver(wd_tag, gm_dict);

    # Users may run the model using prescribed weather drivers by setting p_on, t_on, θ_on to false, or run the energy/water budget by setting those to true.
    df_pres = EmeraldFrontier.simulation!(wd_tag, gm_dict; appending=true, selection = 1:240, p_on = false, t_on = false, θ_on = false);
    df_simu = EmeraldFrontier.simulation!(wd_tag, gm_dict; appending=true, selection = 1:240, p_on = true, t_on = true, θ_on = true);
    @test true;
end
