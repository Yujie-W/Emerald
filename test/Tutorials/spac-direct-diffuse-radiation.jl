using Test

import Emerald.Namespace as ENS
import Emerald.SPAC as ESPAC
import Emerald.Land as ELAND


config = ENS.SPACConfig(Float64);
spac = ENS.BulkSPAC(config);


ESPAC.initialize_spac!(config, spac);
ESPAC.spac!(config, spac, 7200);
ELAND.GPP(config, spac)


spac.meteo.rad_sw.e_dir .*= 0.8;
spac.meteo.rad_sw.e_dif .*= 1.2;
ESPAC.spac!(config, spac, 7200);
ELAND.GPP(config, spac)


spac.meteo.rad_sw.e_dir[600 .<= config.CONSTANTS.SPECTRA.Λ .<= 700] .*= 0.1;
ESPAC.spac!(config, spac, 7200);
ELAND.GPP(config, spac)
