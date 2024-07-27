# This file constains function to write the photosynthetic rates into leaf flux struct

#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-25: add leaf_photosynthesis! function for canopy layer
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(leaf::CanopyLayer{FT}, air::AirLayer{FT}, mode::Union{GCO₂Mode, PCO₂Mode}; rd_only::Bool = false) where {FT}

Updates leaf photosynthetic rates for the leaf based on leaf stomtal model, given
- `leaf` `CanopyLayer` type structure
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `GCO₂Mode` or `PCO₂Mode` to determine whether to use CO₂ partial pressure or concentration to compute photosynthetic rates
- `rd_only` Whether to compute respiration rate only

"""
function leaf_photosynthesis! end;

# This method takes out stomtal model out and use it to determine whether to apply beta to Vcmax, Jmax, and Rd
leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::CanopyLayer{FT},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode};
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(cache, leaf, air, mode, leaf.flux.trait.stomatal_model; rd_only = rd_only);

# if stomtal model is not empirical model, then use the default β = 1
leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::CanopyLayer{FT},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            sm::AbstractStomataModel{FT};
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(cache, leaf, air, mode, FT(1); rd_only = rd_only);

# if stomtal model is empirical model, then determine the β based on the parameter Y (if Vcmax, scale Vcmax, Jmax, and Rd)
leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::CanopyLayer{FT},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}};
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(cache, leaf, air, mode, sm.β, sm.β.PARAM_Y; rd_only = rd_only);

leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::CanopyLayer{FT},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            β::BetaFunction{FT},
            param_y::BetaParameterG1;
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(cache, leaf, air, mode, FT(1); rd_only = rd_only);

leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::CanopyLayer{FT},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            β::BetaFunction{FT},
            param_y::BetaParameterVcmax;
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(cache, leaf, air, mode, leaf.flux.auxil.β; rd_only = rd_only);

# This method computes and save the photosynthetic rates into leaf flux struct for GCO₂Mode
leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::CanopyLayer{FT},
            air::AirLayer{FT},
            mode::GCO₂Mode,
            β::FT;
            rd_only::Bool = false) where {FT} = (
    if rd_only
        leaf.photosystem.auxil.r_d  = leaf.photosystem.trait.r_d25 * temperature_correction(leaf.photosystem.trait.TD_R, leaf.energy.s_aux.t);
        leaf.flux.auxil.a_n        .= -leaf.photosystem.auxil.r_d;
        leaf.flux.auxil.a_g        .= 0;
        leaf.flux.auxil.etr        .= 0;
        leaf.flux.auxil.ϕ_f_shaded  = 0;
        leaf.flux.auxil.ϕ_f_sunlit .= 0;
        leaf.flux.auxil.ϕ_f1       .= 0;
        leaf.flux.auxil.ϕ_f2       .= 0;
        leaf.flux.auxil.ϕ_p        .= 0;

        return nothing
    end;

    photosystem_temperature_dependence!(leaf.photosystem, air, leaf.energy.s_aux.t);
    photosystem_electron_transport!(cache, leaf.photosystem, leaf.flux.auxil.ppar, leaf.flux.auxil.p_CO₂_i; β = β);
    rubisco_limited_rate!(cache, leaf.photosystem, air, leaf.flux.auxil.g_CO₂; β = β);
    light_limited_rate!(cache, leaf.photosystem, air, leaf.flux.auxil.g_CO₂; β = β);
    product_limited_rate!(leaf.photosystem, air, leaf.flux.auxil.g_CO₂; β = β);
    colimit_photosynthesis!(leaf.photosystem; β = β);

    # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
    leaf.flux.auxil.p_CO₂_i .= air.s_aux.ps[2] .- leaf.photosystem.auxil.a_n ./ leaf.flux.auxil.g_CO₂   .* air.state.p_air .* FT(1e-6);
    leaf.flux.auxil.p_CO₂_s .= air.s_aux.ps[2] .- leaf.photosystem.auxil.a_n ./ leaf.flux.auxil.g_CO₂_b .* air.state.p_air .* FT(1e-6);

    # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
    photosystem_electron_transport!(cache, leaf.photosystem, leaf.flux.auxil.ppar, leaf.flux.auxil.p_CO₂_i; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(cache, leaf.photosystem, leaf.flux.auxil.ppar; β = β);

    # save the rates and to leaf (copy here because the photosystem auxil valuse would change when updating stomatal conductance)
    leaf.flux.auxil.a_n .= leaf.photosystem.auxil.a_n;
    leaf.flux.auxil.a_g .= leaf.photosystem.auxil.a_g;
    leaf.flux.auxil.etr .= leaf.photosystem.auxil.a_g ./ leaf.photosystem.auxil.e2c;
    leaf.flux.auxil.ϕ_p .= leaf.photosystem.auxil.ϕ_p;
    leaf.flux.auxil.ϕ_f_shaded = leaf.photosystem.auxil.ϕ_f[end];
    # leaf.flux.auxil.ϕ_f_sunlit[:] .= view(leaf.photosystem.auxil.ϕ_f, 1:length(leaf.photosystem.auxil.ϕ_f)-1);
    nrow = size(leaf.flux.auxil.ϕ_f_sunlit, 1);
    for j in axes(leaf.flux.auxil.ϕ_f_sunlit, 2)
        leaf.flux.auxil.ϕ_f_sunlit[:,j] .= view(leaf.photosystem.auxil.ϕ_f, ((j-1)*nrow+1):(j*nrow));
    end;

    return nothing
);

# This method computes and save the photosynthetic rates into leaf flux struct for PCO₂Mode
leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::CanopyLayer{FT},
            air::AirLayer{FT},
            mode::PCO₂Mode,
            β::FT;
            rd_only::Bool = false) where {FT} = (
    if rd_only
        leaf.photosystem.auxil.r_d   = leaf.photosystem.trait.r_d25 * temperature_correction(leaf.photosystem.trait.TD_R, leaf.energy.s_aux.t);
        leaf.flux.auxil.a_n_sunlit  .= -leaf.photosystem.auxil.r_d;
        leaf.flux.auxil.a_g_sunlit  .= 0;
        leaf.flux.auxil.etr_sunlit  .= 0;
        leaf.flux.auxil.ϕ_f_sunlit  .= 0;
        leaf.flux.auxil.ϕ_f1_sunlit .= 0;
        leaf.flux.auxil.ϕ_f2_sunlit .= 0;
        leaf.flux.auxil.a_n_shaded   = -leaf.photosystem.auxil.r_d;
        leaf.flux.auxil.a_g_shaded   = 0;
        leaf.flux.auxil.etr_shaded   = 0;
        leaf.flux.auxil.ϕ_f_shaded   = 0;
        leaf.flux.auxil.ϕ_f1_shaded  = 0;
        leaf.flux.auxil.ϕ_f2_shaded  = 0;

        return nothing
    end;

    photosystem_temperature_dependence!(leaf.photosystem, air, leaf.energy.s_aux.t);

    # loop through the ppars for sunlit leaf
    for i in eachindex(leaf.flux.auxil.ppar_sunlit)
        # calculate the photosynthetic rates
        photosystem_electron_transport!(leaf.photosystem, leaf.flux.auxil.ppar_sunlit[i], leaf.flux.auxil.p_CO₂_i_sunlit[i]; β = β);
        rubisco_limited_rate!(leaf.photosystem, leaf.flux.auxil.p_CO₂_i_sunlit[i]; β = β);
        light_limited_rate!(leaf.photosystem);
        product_limited_rate!(leaf.photosystem, leaf.flux.auxil.p_CO₂_i_sunlit[i]; β = β);
        colimit_photosynthesis!(leaf.photosystem; β = β);

        # update the fluorescence related parameters
        photosystem_coefficients!(leaf.photosystem, leaf.flux.auxil.ppar_sunlit[i]; β = β);

        # save the rates and to leaf
        leaf.flux.auxil.a_n_sunlit[i]  = leaf.photosystem.auxil.a_n;
        leaf.flux.auxil.a_g_sunlit[i]  = leaf.photosystem.auxil.a_g;
        leaf.flux.auxil.etr_sunlit[i]  = leaf.photosystem.auxil.a_g / leaf.photosystem.auxil.e2c;
        leaf.flux.auxil.ϕ_d_sunlit[i]  = leaf.photosystem.auxil.ϕ_d;
        leaf.flux.auxil.ϕ_f_sunlit[i]  = leaf.photosystem.auxil.ϕ_f;
        leaf.flux.auxil.ϕ_n_sunlit[i]  = leaf.photosystem.auxil.ϕ_n;
        leaf.flux.auxil.ϕ_p_sunlit[i]  = leaf.photosystem.auxil.ϕ_p;
        leaf.flux.auxil.ϕ_f1_sunlit[i] = leaf.photosystem.auxil.ϕ_f1;
        leaf.flux.auxil.ϕ_f2_sunlit[i] = leaf.photosystem.auxil.ϕ_f2;
    end;

    # run the model for shaded leaf
    photosystem_electron_transport!(leaf.photosystem, leaf.flux.auxil.ppar_shaded, leaf.flux.auxil.p_CO₂_i_shaded; β = β);
    rubisco_limited_rate!(leaf.photosystem, leaf.flux.auxil.p_CO₂_i_shaded; β = β);
    light_limited_rate!(leaf.photosystem);
    product_limited_rate!(leaf.photosystem, leaf.flux.auxil.p_CO₂_i_shaded; β = β);
    colimit_photosynthesis!(leaf.photosystem; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(leaf.photosystem, leaf.flux.auxil.ppar_shaded; β = β);

    # save the rates and to leaf
    leaf.flux.auxil.a_n_shaded = leaf.photosystem.auxil.a_n;
    leaf.flux.auxil.a_g_shaded = leaf.photosystem.auxil.a_g;
    leaf.flux.auxil.etr_shaded = leaf.photosystem.auxil.a_g / leaf.photosystem.auxil.e2c;
    leaf.flux.auxil.ϕ_d_shaded = leaf.photosystem.auxil.ϕ_d;
    leaf.flux.auxil.ϕ_f_shaded = leaf.photosystem.auxil.ϕ_f;
    leaf.flux.auxil.ϕ_n_shaded = leaf.photosystem.auxil.ϕ_n;
    leaf.flux.auxil.ϕ_p_shaded = leaf.photosystem.auxil.ϕ_p;

    return nothing
);
