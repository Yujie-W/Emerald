# This script is meant to speed up the SPAC model

using Emerald
using Revise

FT = Float64;

config = EmeraldLand.Namespace.SPACConfiguration(FT);
spac_c = EmeraldLand.Namespace.BulkSPAC(config; use_leaf = false);
EmeraldLand.SPAC.initialize_spac!(config, spac_c);
EmeraldLand.SPAC.spac!(config, spac_c, FT(1));

@time EmeraldLand.SPAC.spac!(config, spac_c, FT(3600));

@info "GPP and SIF" EmeraldLand.SPAC.GPP(spac_c) EmeraldLand.SPAC.TROPOMI_SIF740(config, spac_c);


spac = EmeraldLand.Namespace.BulkSPAC(config; use_leaf = false);
@time EmeraldLand.SPAC.initialize_spac!(config, spac);  # allocations because of setfield! for numbers
@time EmeraldLand.SPAC.spac!(config, spac, FT(3600));   # 1 allocation of 16 bytes (in the adjust_time function), ignored here
