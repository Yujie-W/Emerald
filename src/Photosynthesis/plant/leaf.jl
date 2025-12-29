# This method takes out stomtal model out and use it to determine whether to apply beta to Vcmax, Jmax, and Rd
leaf_photosynthesis!(
            config::SPACConfig{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT};
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(config, leaf, air, config.METHODS.STOMATAL_MODEL; rd_only = rd_only);

# if stomtal model is not empirical model, then use the default β = 1
leaf_photosynthesis!(
            config::SPACConfig{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT},
            sm::AbstractStomatalConductanceModel{FT};
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(config, leaf, air, FT(1); rd_only = rd_only);

# if stomtal model is empirical model, then determine the β based on the parameter Y (if Vcmax, scale Vcmax, Jmax, and Rd)
leaf_photosynthesis!(
            config::SPACConfig{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT},
            sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}};
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(config, leaf, air, sm.β, sm.β.PARAM_Y; rd_only = rd_only);

leaf_photosynthesis!(
            config::SPACConfig{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT},
            β::BetaFunction{FT},
            param_y::BetaParameterG1;
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(config, leaf, air, FT(1); rd_only = rd_only);

leaf_photosynthesis!(
            config::SPACConfig{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT},
            β::BetaFunction{FT},
            param_y::BetaParameterVcmax;
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(config, leaf, air, leaf.flux.auxil.β; rd_only = rd_only);

# This method computes and save the photosynthetic rates into leaf flux struct for Conductance mode
leaf_photosynthesis!(
            config::SPACConfig{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT},
            β::FT;
            rd_only::Bool = false) where {FT} = (
    if rd_only
        return photosystem_temperature_dependence!(config, leaf.photosystem, air, leaf.energy.auxil.t);
    end;

    photosystem_temperature_dependence!(config, leaf.photosystem, air, leaf.energy.auxil.t);
    photosystem_electron_transport!(config, leaf.photosystem, leaf.flux.auxil.ppar, leaf.flux.auxil.p_CO₂_i; β = β);
    rubisco_limited_rate!(config, leaf.photosystem, air, leaf.flux.auxil.g_CO₂; β = β);
    light_limited_rate!(config, leaf.photosystem, air, leaf.flux.auxil.g_CO₂; β = β);
    product_limited_rate!(config, leaf.photosystem, air, leaf.flux.auxil.g_CO₂; β = β);
    colimit_photosynthesis!(config, leaf.photosystem; β = β);

    # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
    leaf.flux.auxil.p_CO₂_i = air.auxil.ps[2] - leaf.photosystem.auxil.a_n / leaf.flux.auxil.g_CO₂   * air.state.p_air * FT(1e-6);
    leaf.flux.auxil.p_CO₂_s = air.auxil.ps[2] - leaf.photosystem.auxil.a_n / leaf.flux.auxil.g_CO₂_b * air.state.p_air * FT(1e-6);

    # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
    photosystem_electron_transport!(config, leaf.photosystem, leaf.flux.auxil.ppar, leaf.flux.auxil.p_CO₂_i; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(config, leaf.photosystem, leaf.flux.auxil.ppar; β = β);

    return nothing
);
