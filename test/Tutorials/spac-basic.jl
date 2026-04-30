using Test

import Emerald.Namespace as ENS
import Emerald.SPAC as ESPAC
import Emerald.Land as ELAND


config = ENS.SPACConfig(Float64);
spac = ENS.BulkSPAC(config);
ESPAC.initialize_spac!(config, spac);
spac.plant.pool.c_pool = Inf;
ESPAC.prescribe_traits!(config, spac; sai = 0, lai = 3);
ESPAC.spac!(config, spac, 0);


for lai in 0.1:0.1:6.0
    ESPAC.prescribe_traits!(config, spac; sai = 0, lai = lai);
    ESPAC.spac!(config, spac, 0);
    @show lai ELAND.SHORTWAVE_OUT(config, spac) / ELAND.SHORTWAVE_IN(config, spac);
end;


ESPAC.prescribe_traits!(config, spac; sai = 0, lai = 2);
for sza in 0:5:85
    spac.canopy.sun_geometry.state.sza = sza;
    ESPAC.spac!(config, spac, 0);
    @show sza ELAND.SHORTWAVE_OUT(config, spac) / ELAND.SHORTWAVE_IN(config, spac);
end;
