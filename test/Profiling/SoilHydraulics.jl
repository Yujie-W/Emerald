# make sure the LeafOptics module is fast enough that it won't be a bottleneck in simulations
#
# julia --project --track-allocation=user
#


using Profile

import Emerald.Namespace as ENS
import Emerald.SoilHydraulics as ESH
import Emerald.SPAC as ESPAC


config = ENS.SPACConfig(Float64);
spac = ENS.BulkSPAC(config);
ESPAC.initialize_spac!(config, spac);
ESPAC.spac!(config, spac, 1);

@time ESH.soil_profiles!(config, spac);
@time ESH.soil_budgets!(config, spac, 1.0);

Profile.clear_malloc_data();

@time ESH.soil_profiles!(config, spac);
@time ESH.soil_budgets!(config, spac, 1.0);


exit()
