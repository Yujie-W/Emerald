# Test the Leaf and CanopyLayer photosynthesis in parallel
using Emerald


# define the SPACs
FT = Float64;
config = EmeraldLand.Namespace.SPACConfiguration(FT);
spac_c = EmeraldLand.Namespace.BulkSPAC(config);
EmeraldLand.SPAC.initialize_spac!(config, spac_c);


# get the leaves and air
leaf_c = spac_c.plant.leaves[1];
air = spac_c.airs[1];

EmeraldLand.Photosynthesis.photosystem_temperature_dependence!(leaf_c.photosystem, air, leaf_c.energy.s_aux.t);
