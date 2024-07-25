# This script is meant to speed up the SPAC model

using Emerald
using Revise

FT = Float64;

config = EmeraldLand.Namespace.SPACConfiguration(FT);

spac = EmeraldLand.Namespace.BulkSPAC(config);
EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, FT(1));


spac = EmeraldLand.Namespace.BulkSPAC(config);

@time EmeraldLand.SPAC.initialize_spac!(config, spac);  # allocations because of setfield! for numbers

@time EmeraldLand.SPAC.spac!(config, spac, FT(3600));   # 1 allocation of 16 bytes (in the adjust_time function), ignored here
