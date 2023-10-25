using Emerald;


wdtag = "wd1";
gmdict = EmeraldFrontier.gm_dict(EmeraldFrontier.GriddingMachineLabels(year=2019), -2.6, -52.9);
EmeraldFrontier.simulation!(wdtag, gmdict; appending=false, selection = 1:480, saving = "/home/wyujie/ProjectData/2023_Cytochrome_Model/test-run-c3vjp.nc");
