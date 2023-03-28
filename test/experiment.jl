#
# Experiment of whether to prescribe water and energy related parameters
#
using Emerald;

wdtag = "wd1";
gmdict = EmeraldFrontier.gm_dict(EmeraldFrontier.GriddingMachineLabels(year=2019), 35.1, 115.2);
EmeraldFrontier.simulation!(wdtag, gmdict; appending=true, p_on = false, t_on = false, θ_on = false, saving = "simulation.priscription.nc");
EmeraldFrontier.simulation!(wdtag, gmdict; appending=true, p_on = true, t_on = true, θ_on = true, saving = "simulation.simulation.nc");
