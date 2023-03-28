#
# Experiment of whether to prescribe water and energy related parameters
#
using Emerald;

wd_tag = "wd1";
gm_dict = EmeraldFrontier.gm_dict(EmeraldFrontier.GriddingMachineLabels(year=2019), 35.1, 115.2);
df_pres = EmeraldFrontier.simulation!(wd_tag, gm_dict; appending=true, selection = 1:240, p_on = false, t_on = false, θ_on = false);
df_simu = EmeraldFrontier.simulation!(wd_tag, gm_dict; appending=true, selection = 1:240, p_on = true, t_on = true, θ_on = true);
