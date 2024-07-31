import Emerald.EmeraldData.GlobalDatasets as GD
import Emerald.EmeraldFrontier as EF


gm_tag = "gm2";
wd_tag = "wd1";




# the case using DataFrame
@time gm_dict = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2019), 38.74, -92.20);
@time df_simu = EF.simulation!(wd_tag, gm_dict; appending = false, selection = 1:24);




# the case using NamedTuple
@time gm_dict = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2019), 38.74, -92.20);
@time df_simu = EF.simulation_nt!(wd_tag, gm_dict; appending = false, selection = 1:24);
