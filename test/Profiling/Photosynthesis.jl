# make sure the LeafOptics module is fast enough that it won't be a bottleneck in simulations
#
# julia --project --track-allocation=user
#
using BenchmarkTools
using Profile

import Emerald.Namespace as ENS
import Emerald.Photosynthesis as EPH
import Emerald.SPAC as ESPAC


config = ENS.SPACConfig(Float64);
spac = ENS.BulkSPAC(config);
ESPAC.initialize_spac!(config, spac);
lbio = spac.plant.leaves[end].bio;

@time EPH.plant_photosynthesis!(config, spac);

Profile.clear_malloc_data();

@time EPH.plant_photosynthesis!(config, spac);
