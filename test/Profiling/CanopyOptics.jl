# make sure the LeafOptics module is fast enough that it won't be a bottleneck in simulations
#
# julia --project --track-allocation=user
#


using Profile

import Emerald.Namespace as ENS
import Emerald.CanopyOptics as ECO
import Emerald.SPAC as ESPAC


config = ENS.SPACConfig(Float64);
spac = ENS.BulkSPAC(config);
ESPAC.initialize_spac!(config, spac);

@time ECO.canopy_radiation!(config, spac);

Profile.clear_malloc_data();

@time ECO.canopy_radiation!(config, spac);


exit()
