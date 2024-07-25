# This script is meant to speed up the SPAC model

using Emerald
using Revise

FT = Float64;

config = EmeraldLand.Namespace.SPACConfiguration(FT);


spac_l = EmeraldLand.Namespace.BulkSPAC(config; use_leaf = true);
spac_c = EmeraldLand.Namespace.BulkSPAC(config; use_leaf = false);
EmeraldLand.SPAC.initialize_spac!(config, spac_l);
EmeraldLand.SPAC.initialize_spac!(config, spac_c);

@time EmeraldLand.SPAC.spac!(config, spac_l, FT(3600));
@time EmeraldLand.SPAC.spac!(config, spac_c, FT(3600));

@info "GPP" EmeraldLand.SPAC.GPP(spac_c) EmeraldLand.SPAC.GPP(spac_l);


spac = EmeraldLand.Namespace.BulkSPAC(config; use_leaf = true);
@time EmeraldLand.SPAC.initialize_spac!(config, spac);  # allocations because of setfield! for numbers
@time EmeraldLand.SPAC.spac!(config, spac, FT(3600));   # 1 allocation of 16 bytes (in the adjust_time function), ignored here
