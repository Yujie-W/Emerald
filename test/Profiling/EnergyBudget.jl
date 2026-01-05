# make sure the LeafOptics module is fast enough that it won't be a bottleneck in simulations
#
# julia --project --track-allocation=user
#


using BenchmarkTools
using Profile

import Emerald.Namespace as ENS
import Emerald.EnergyBudget as EEB
import Emerald.SPAC as ESPAC


config = ENS.SPACConfig(Float64);
spac = ENS.BulkSPAC(config);
ESPAC.initialize_spac!(config, spac);
ESPAC.spac!(config, spac, 1);

@time EEB.spac_energy_flow!(config, spac);
@time EEB.spac_energy_budget!(spac, 1.0);

Profile.clear_malloc_data();

@time EEB.spac_energy_flow!(config, spac);
@time EEB.spac_energy_budget!(spac, 1.0);


exit()
