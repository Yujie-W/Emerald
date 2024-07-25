# Test the Leaf and CanopyLayer photosynthesis in parallel
using Emerald


# define the SPACs
FT = Float64;
config = EmeraldLand.Namespace.SPACConfiguration(FT);
spac_c = EmeraldLand.Namespace.BulkSPAC(config; use_leaf = false);
spac_l = EmeraldLand.Namespace.BulkSPAC(config; use_leaf = true);
EmeraldLand.SPAC.initialize_spac!(config, spac_c);
EmeraldLand.SPAC.initialize_spac!(config, spac_l);


# get the leaves and air
leaf_c = spac_c.plant.leaves[1];
leaf_l = spac_l.plant.leaves[1];
air = spac_c.airs[1];

EmeraldLand.Photosynthesis.photosystem_temperature_dependence!(leaf_c.photosystem, air, leaf_c.energy.s_aux.t);
EmeraldLand.Photosynthesis.photosystem_temperature_dependence!(leaf_l.photosystem, air, leaf_l.energy.s_aux.t);


# Compare the electron transport
ppars = ones(FT, length(leaf_c.photosystem.auxil.a_n)) .* 1000;
pis = ones(FT, length(leaf_c.photosystem.auxil.a_n)) .* 20;
glcs = ones(FT, length(leaf_c.photosystem.auxil.a_n)) .* 0.1;
EmeraldLand.Photosynthesis.photosystem_electron_transport!(leaf_c.photosystem, ppars, pis);
EmeraldLand.Photosynthesis.photosystem_electron_transport!(leaf_l.photosystem, ppars[1], pis[1]);
@info "Debugging" leaf_c.photosystem.auxil.e2c leaf_l.photosystem.auxil.e2c leaf_c.photosystem.auxil.j leaf_l.photosystem.auxil.j;


# Compare the Ac
EmeraldLand.Photosynthesis.rubisco_limited_rate!(leaf_c.photosystem, pis);
EmeraldLand.Photosynthesis.rubisco_limited_rate!(leaf_l.photosystem, pis[1]);
@info "Debugging" leaf_c.photosystem.auxil.a_c leaf_l.photosystem.auxil.a_c;


# Compare the Aj
EmeraldLand.Photosynthesis.light_limited_rate!(leaf_c.photosystem);
EmeraldLand.Photosynthesis.light_limited_rate!(leaf_l.photosystem);
@info "Debugging" leaf_c.photosystem.auxil.a_j leaf_l.photosystem.auxil.a_j;


# Compare the Ap
EmeraldLand.Photosynthesis.product_limited_rate!(leaf_c.photosystem, pis);
EmeraldLand.Photosynthesis.product_limited_rate!(leaf_l.photosystem, pis[1]);
@info "Debugging" leaf_c.photosystem.auxil.a_p leaf_l.photosystem.auxil.a_p;


# Compare the Ag and An
EmeraldLand.Photosynthesis.colimit_photosynthesis!(leaf_c.photosystem);
EmeraldLand.Photosynthesis.colimit_photosynthesis!(leaf_l.photosystem);
@info "Debugging" leaf_c.photosystem.auxil.a_g leaf_l.photosystem.auxil.a_g;
@info "Debugging" leaf_c.photosystem.auxil.a_n leaf_l.photosystem.auxil.a_n;


# Compare the model at G Mode
EmeraldLand.Photosynthesis.rubisco_limited_rate!(leaf_c.photosystem, air, glcs);
EmeraldLand.Photosynthesis.rubisco_limited_rate!(leaf_l.photosystem, air, glcs[1]);
@info "Debugging" leaf_c.photosystem.auxil.a_c leaf_l.photosystem.auxil.a_c;

EmeraldLand.Photosynthesis.light_limited_rate!(leaf_c.photosystem, air, glcs);
EmeraldLand.Photosynthesis.light_limited_rate!(leaf_l.photosystem, air, glcs[1]);
@info "Debugging" leaf_c.photosystem.auxil.a_j leaf_l.photosystem.auxil.a_j;

EmeraldLand.Photosynthesis.product_limited_rate!(leaf_c.photosystem, air, glcs);
EmeraldLand.Photosynthesis.product_limited_rate!(leaf_l.photosystem, air, glcs[1]);
@info "Debugging" leaf_c.photosystem.auxil.a_p leaf_l.photosystem.auxil.a_p;

EmeraldLand.Photosynthesis.colimit_photosynthesis!(leaf_c.photosystem);
EmeraldLand.Photosynthesis.colimit_photosynthesis!(leaf_l.photosystem);
@info "Debugging" leaf_c.photosystem.auxil.a_g leaf_l.photosystem.auxil.a_g;
@info "Debugging" leaf_c.photosystem.auxil.a_n leaf_l.photosystem.auxil.a_n;
