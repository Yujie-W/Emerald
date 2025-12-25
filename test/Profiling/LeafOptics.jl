# make sure the LeafOptics module is fast enough that it won't be a bottleneck in simulations
#
# julia --project --track-allocation=user
#
using BenchmarkTools
using Profile

import Emerald.LeafOptics as ELO
import Emerald.Namespace as ENS
import Emerald.SPAC as ESPAC


config = ENS.SPACConfig(Float64);
spac = ENS.BulkSPAC(config);
ESPAC.initialize_spac!(config, spac);
lbio = spac.plant.leaves[end].bio;

@time ELO.leaf_spectra!(config, spac.plant.leaves[end].bio, spac.cache, 5.0);
@time ELO.plant_leaf_spectra!(config, spac);

Profile.clear_malloc_data()

@time ELO.leaf_spectra!(config, lbio, spac.cache, 5.0);
@time ELO.plant_leaf_spectra!(config, spac);
