#
# Experiment to test global simulation
#
using Emerald;

FT = Float64;
@time EmeraldEarth.add_threads!(20, FT);

@time dts = EmeraldEarth.LandDatasets{FT}("gm2", 2020);
@time mat = EmeraldEarth.gm_grids(dts);
@time sts = Matrix{Union{Nothing,EmeraldLand.Namespace.MultiLayerSPACState{FT}}}(nothing, size(dts.t_lm));
@time wdr = EmeraldData.ERA5.ERA5SingleLevelsDriver();

# debug the code on a single core
#=
@time sts = Matrix{Union{Nothing,EmeraldLand.Namespace.MultiLayerSPACState{FT}}}(nothing, size(dts.t_lm));
@time wd1 = EmeraldEarth.wd_grids(dts, wdr, 1);
@time sts = EmeraldEarth.simulation!(mat, wd1, sts);
=#

# for debugging use
rm("test.nc"; force = true);
for ind in 1:6:24
    @time wdx = EmeraldEarth.wd_grids(dts, wdr, ind);
    @time sts = EmeraldEarth.simulation!(mat, wdx, sts);
    @time EmeraldEarth.save_simulations!("test.nc", sts, ind);
end;

# visualize the simulations to make sure the simulation is okay
EmeraldVisualization.animate_nc!("test.nc", "BETA"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "CSIF"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "ETR"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "GPP"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "EVI"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "NDVI"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "NIRv"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "PAR"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "PPAR"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "SIF₆₈₃"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "SIF₇₄₀"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "SIF₇₅₉"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "SIF₇₇₀"; filename = "test.gif", fps = 2);
EmeraldVisualization.animate_nc!("test.nc", "TRANSPIRATION"; filename = "test.gif", fps = 2);
