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

    return nothing
);
