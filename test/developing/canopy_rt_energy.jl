# whether energy is conserved

using Emerald


FT = Float64;

config = EmeraldLand.Namespace.SPACConfiguration(FT);
spac = EmeraldLand.Namespace.BulkSPAC(config);
EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, FT(3600));


# incoming radiation
e_incoming = (spac.meteo.rad_sw.e_dir + spac.meteo.rad_sw.e_dif)' * config.SPECTRA.ΔΛ / 1000;

# total absorption
e_absorbed = sum(spac.canopy.sun_geometry.auxil.r_net_sw_leaf .+ spac.canopy.sun_geometry.auxil.r_net_sw_stem) + spac.soil_bulk.auxil.r_net_sw;

# outgoing radiation
e_outgoing = spac.canopy.sun_geometry.auxil.e_difꜛ[:,1]' * config.SPECTRA.ΔΛ / 1000;

@info "Energy balance (W m⁻²)" e_incoming e_absorbed e_outgoing e_incoming - e_absorbed - e_outgoing;
