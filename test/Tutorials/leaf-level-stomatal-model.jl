import Emerald.Namespace as ENS
import Emerald.ResearchTools as ERT
import Emerald.SPAC as ESPAC


config_c3 = ERT.LeafLevelSetup.leaf_level_config(Float64);
config_c4 = ERT.LeafLevelSetup.leaf_level_config(Float64);
config_c4_new = ERT.LeafLevelSetup.leaf_level_config(Float64);
config_c4_new.FEATURES.NEW_C4_STOMATAL_METHODS = true;

c3leaf = ERT.LeafLevelSetup.leaf_level_leaf(config_c3, "C3");
c4leaf = ERT.LeafLevelSetup.leaf_level_leaf(config_c4, "C4");
c4leaf_new = ERT.LeafLevelSetup.leaf_level_leaf(config_c4_new, "C4");
air = ENS.AirLayer{Float64}();

c3leaf.flux.auxil.ppar .= 1500;
c3leaf.energy.auxil.t = 305.15;
c3leaf.xylem.auxil.pressure[1] = -2.0;
c4leaf.flux.auxil.ppar .= 1500;
c4leaf.energy.auxil.t = 305.15;
c4leaf.xylem.auxil.pressure[1] = -2.0;
c4leaf_new.flux.auxil.ppar .= 1500;
c4leaf_new.energy.auxil.t = 305.15;
c4leaf_new.xylem.auxil.pressure[1] = -2.0;
ESPAC.prescribe_air!(air; f_CO₂ = 400, vpd = 3000);


c3leaf.flux.state.g_H₂O_s .= 0.12;
c4leaf.flux.state.g_H₂O_s .= 0.12;
c4leaf_new.flux.state.g_H₂O_s .= 0.12;
ERT.Stomata.∂A∂E_∂Θ∂E(config_c3, c3leaf, air)
ERT.Stomata.∂A∂E_∂Θ∂E(config_c4, c4leaf, air)
ERT.Stomata.∂A∂E_∂Θ∂E(config_c4_new, c4leaf_new, air)


ERT.Stomata.steady_state_gs!(config_c3, c3leaf, air);
c3leaf.flux.state.g_H₂O_s[1]
c3leaf.flux.auxil.∂A∂E[1]
c3leaf.flux.auxil.∂Θ∂E[1]
ERT.Stomata.steady_state_gs!(config_c4, c4leaf, air);
c4leaf.flux.state.g_H₂O_s[1]
c4leaf.flux.auxil.∂A∂E[1]
c4leaf.flux.auxil.∂Θ∂E[1]
ERT.Stomata.steady_state_gs!(config_c4_new, c4leaf_new, air);
c4leaf_new.flux.state.g_H₂O_s[1]
c4leaf_new.flux.auxil.∂A∂E[1]
c4leaf_new.flux.auxil.∂Θ∂E[1]
