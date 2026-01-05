# make sure the LeafOptics module is fast enough that it won't be a bottleneck in simulations
#
# julia --project --track-allocation=user
#


using Profile

import Emerald.Namespace as ENS
import Emerald.Photosynthesis as EPH
import Emerald.SPAC as ESPAC


config = ENS.SPACConfig(Float64);
cyto_model = false;
if cyto_model
    config.METHODS.C3_AC_METHOD = ENS.AcMethodC3VcmaxPi();
    config.METHODS.C3_AJ_METHOD = ENS.AjMethodC3VqmaxPi();
    config.METHODS.C3_AP_METHOD = ENS.ApMethodC3Vcmax();
    config.METHODS.COLIMIT_J = ENS.SerialColimit();
    config.METHODS.FLUORESCENCE_METHOD = ENS.CytochromeFluorescenceModel();
end;
spac = ENS.BulkSPAC(config);
ESPAC.initialize_spac!(config, spac);

@time EPH.plant_photosynthesis!(config, spac);
@time EPH.plant_carbon_budget!(spac, 1.0);

Profile.clear_malloc_data();

@time EPH.plant_photosynthesis!(config, spac);
@time EPH.plant_carbon_budget!(spac, 1.0);


exit()
