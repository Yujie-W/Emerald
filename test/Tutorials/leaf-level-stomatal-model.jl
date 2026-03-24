import Emerald.Namespace as ENS
import Emerald.ResearchTools as ERT
import Emerald.SPAC as ESPAC




config_c4_new = ERT.LeafLevelSetup.leaf_level_config(Float64);
config_c4_new.FEATURES.NEW_C4_STOMATAL_METHODS = true;
config_c4_new.METHODS.C4_AP_METHOD = ENS.ApMethodC4VpmaxPi();
config_c4_new.METHODS.STOMATAL_MODEL = ENS.Sperry2SM{Float64}();


# test the VPD responses
for vpd in collect(200:100:3000)
    c4leaf_new = ERT.LeafLevelSetup.leaf_level_leaf(config_c4_new, "C4");
    air = ENS.AirLayer{Float64}();
    c4leaf_new.flux.auxil.ppar .= 1500;
    c4leaf_new.energy.auxil.t = 300.15;
    c4leaf_new.xylem.auxil.pressure[1] = -2.0;
    ESPAC.prescribe_air!(air; f_CO₂ = 400, vpd = vpd);
    ERT.Stomata.steady_state_gs!(config_c4_new, c4leaf_new, air);
    am = min.(c4leaf_new.photosystem.auxil.a_p, c4leaf_new.photosystem.auxil.a_j) .- c4leaf_new.photosystem.auxil.r_d;
    @info "VPD = $vpd" c4leaf_new.flux.state.g_H₂O_s[1] c4leaf_new.flux.auxil.p_CO₂_i[1] am[1];
end;





begin
    c4leaf_new = ERT.LeafLevelSetup.leaf_level_leaf(config_c4_new, "C4");
    air = ENS.AirLayer{Float64}();
    c4leaf_new.flux.auxil.ppar .= 1500;
    c4leaf_new.energy.auxil.t = 300.15;
    c4leaf_new.xylem.auxil.pressure[1] = -2.0;
    ESPAC.prescribe_air!(air; f_CO₂ = 400, vpd = 2400);
    ERT.Stomata.steady_state_gs!(config_c4_new, c4leaf_new, air);
    am = min.(c4leaf_new.photosystem.auxil.a_p, c4leaf_new.photosystem.auxil.a_j) .- c4leaf_new.photosystem.auxil.r_d;
    @info "VPD = 2400" c4leaf_new.flux.state.g_H₂O_s[1] c4leaf_new.flux.auxil.p_CO₂_i[1] am[1];
end;



# config_c3 = ERT.LeafLevelSetup.leaf_level_config(Float64);
# config_c4 = ERT.LeafLevelSetup.leaf_level_config(Float64);
# config_c4_new = ERT.LeafLevelSetup.leaf_level_config(Float64);
# config_c4_new.FEATURES.NEW_C4_STOMATAL_METHODS = true;


# c3leaf = ERT.LeafLevelSetup.leaf_level_leaf(config_c3, "C3");
# c4leaf = ERT.LeafLevelSetup.leaf_level_leaf(config_c4, "C4");

# c3leaf.flux.auxil.ppar .= 1500;
# c3leaf.energy.auxil.t = 305.15;
# c3leaf.xylem.auxil.pressure[1] = -2.0;
# c4leaf.flux.auxil.ppar .= 1500;
# c4leaf.energy.auxil.t = 305.15;
# c4leaf.xylem.auxil.pressure[1] = -2.0;


# c3leaf.flux.state.g_H₂O_s .= 0.12;
# c4leaf.flux.state.g_H₂O_s .= 0.12;
# c4leaf_new.flux.state.g_H₂O_s .= 0.12;
# ERT.Stomata.∂A∂E_∂Θ∂E(config_c3, c3leaf, air)
# ERT.Stomata.∂A∂E_∂Θ∂E(config_c4, c4leaf, air)
# ERT.Stomata.∂A∂E_∂Θ∂E(config_c4_new, c4leaf_new, air)


# ERT.Stomata.steady_state_gs!(config_c3, c3leaf, air);
# c3leaf.flux.state.g_H₂O_s[1]
# c3leaf.flux.auxil.∂A∂E[1]
# c3leaf.flux.auxil.∂Θ∂E[1]
# ERT.Stomata.steady_state_gs!(config_c4, c4leaf, air);
# c4leaf.flux.state.g_H₂O_s[1]
# c4leaf.flux.auxil.∂A∂E[1]
# c4leaf.flux.auxil.∂Θ∂E[1]
# c4leaf_new.flux.state.g_H₂O_s[1]
# c4leaf_new.flux.auxil.∂A∂E[1]
# c4leaf_new.flux.auxil.∂Θ∂E[1]
