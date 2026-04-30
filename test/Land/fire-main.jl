# 1. read in the lat/lon/ratio data
# df = read_csv();

params = [];
for dfr in eachrow(df)
    push!(params, [settings, dfr.lat, dfr.lon, dfr.year, dfr.ratio, filename]);
end;

#=
for p in params
    simulation_yx!(settings, p[1], p[2], p[3]; c3c4 = "C3", saving = p[4]);
end;
=#

using Distributed: pmap, @everywhere
using PkgUtility.DistributedTools: dynamic_workers!


dynamic_workers!(40);
@everywhere include("fire-thread.jl");
pmap(thread_func, params);
