# make sure the LeafOptics module is fast enough that it won't be a bottleneck in simulations
#
# julia --project --track-allocation=user
#


using Profile

import Emerald.Namespace as ENS
import Emerald.SPAC as ESPAC


config = ENS.SPACConfig(Float64);
spac = ENS.BulkSPAC(config);

@time ESPAC.initialize_spac!(config, spac);
@time ESPAC.spac!(config, spac, 1);

Profile.clear_malloc_data();

@time ESPAC.initialize_spac!(config, spac);
@time ESPAC.spac!(config, spac, 1);


exit()
