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

# preloaded driver mode
#=
@time wds = EmeraldEarth.weather_drivers(dts, wdr);
@time wdx = EmeraldEarth.wd_grids(dts, wds, 1);
=#

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




# visualize the simulations to make sure the simulated area
nansts = zeros(Bool, size(sts));
for i in eachindex(sts)
    if !isnothing(sts[i])
        nansts[i] = true;
    end;
end;
@show sum(nansts);
EmeraldVisualization.animate_data!(collect(-179.5:1:180), collect(-89.5:1:90), nansts; filename = "test.png");




# debugging the code using NaNs within the simulations
gpps = EmeraldIO.Netcdf.read_nc("/home/wyujie/Github/Julia/Emerald/test.nc", "NDVI");
nans = [];
for i in 1:360, j in 1:180
    if 1 <= sum(isnan.(gpps[i,j,:])) < size(gpps,3)
        push!(nans, [i,j])
    end;
end;
nans




# for debugging use
@time sts = Matrix{Union{Nothing,EmeraldLand.Namespace.MultiLayerSPACState{FT}}}(nothing, size(dts.t_lm));

@time wd1 = EmeraldEarth.wd_grids(dts, wdr, 1);
@time wd2 = EmeraldEarth.wd_grids(dts, wdr, 4);
@time wd3 = EmeraldEarth.wd_grids(dts, wdr, 7);
@time wd4 = EmeraldEarth.wd_grids(dts, wdr, 10);
@time wd5 = EmeraldEarth.wd_grids(dts, wdr, 13);
@time wd6 = EmeraldEarth.wd_grids(dts, wdr, 16);

iii = 244;
jjj = 120;
sta = sts[iii,jjj];
sta = EmeraldEarth.simulation!(mat[iii,jjj], wd1[iii,jjj], sta);
sta = EmeraldEarth.simulation!(mat[iii,jjj], wd2[iii,jjj], sta);
sta = EmeraldEarth.simulation!(mat[iii,jjj], wd3[iii,jjj], sta);
sta = EmeraldEarth.simulation!(mat[iii,jjj], wd4[iii,jjj], sta);
sta = EmeraldEarth.simulation!(mat[iii,jjj], wd5[iii,jjj], sta);
sta = EmeraldEarth.simulation!(mat[iii,jjj], wd6[iii,jjj], sta);
