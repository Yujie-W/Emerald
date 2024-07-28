import Emerald.EmeraldData.GlobalDatasets as GD
import Emerald.EmeraldFrontier as EF


gm_tag = "gm4";
wd_tag = "wd1";

@time gm_dict = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2019), 38.74, -92.20);

@time df_simu = EF.simulation!(wd_tag, gm_dict; appending = false, selection = 1:24);
