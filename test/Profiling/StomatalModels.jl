# make sure the LeafOptics module is fast enough that it won't be a bottleneck in simulations
#
# julia --project --track-allocation=user
#


using Profile

import Emerald.Namespace as ENS
import Emerald.StomatalModels as ESM
import Emerald.SPAC as ESPAC


config = ENS.SPACConfig(Float64);
spac = ENS.BulkSPAC(config);
ESPAC.initialize_spac!(config, spac);
ESPAC.spac!(config, spac, 1);

@time ESM.stomatal_conductance_profile!(config, spac);
@time ESM.stomatal_conductance!(spac, 1.0);

Profile.clear_malloc_data();

@time ESM.stomatal_conductance_profile!(config, spac);
@time ESM.stomatal_conductance!(spac, 1.0);


exit()
