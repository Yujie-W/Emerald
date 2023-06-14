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

# for debugging use
@time wdx = EmeraldEarth.wd_grids(dts, wdr, 1);
@time sts = EmeraldEarth.simulation!(mat, wdx, sts);
@time EmeraldEarth.save_simulations!("test.nc", sts, 1);
@time wdx = EmeraldEarth.wd_grids(dts, wdr, 2);
@time sts = EmeraldEarth.simulation!(mat, wdx, sts);
@time EmeraldEarth.save_simulations!("test.nc", sts, 2);

# visualize the simulations to make sure the simulated area
nansts = zeros(Bool, size(sts));
for i in eachindex(sts)
    if !isnothing(sts[i])
        nansts[i] = true;
    end;
end;
@show sum(nansts);
EmeraldVisualization.animate_data!(collect(-179.5:1:180), collect(-89.5:1:90), nansts; filename = "test.png");

# function to run the simulation at site level at one time step (parallelizing module)
EmeraldEarth.simulation!(mat[339,36], wdx[339,36])
EmeraldEarth.simulation!(mat[296,120], wdx[296,120])
