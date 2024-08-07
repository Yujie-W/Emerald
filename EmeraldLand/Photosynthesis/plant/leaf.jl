# This file constains function to write the photosynthetic rates into leaf flux struct

#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: set a default p_i from leaf to combine two methods
#     2022-Jan-14: do not update temperature to avoid its impact on plant hydraulics
#     2022-Jan-14: use colimit function to compute a_gross and a_net
#     2022-Jan-14: set a default g_lc from leaf to combine two methods
#     2022-Jan-18: add p_i to electron transport function input variables
#     2022-Feb-28: add support to C3CytochromeModel
#     2022-Jun-27: remove ppar from input variable list of light_limited_rate!
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2022-Jul-07: save a_net and a_gross to Leaf (as PSM may be used for temporary calculations)
#     2022-Jul-12: use β as a must have option (and thus this function becomes a core function of the one above)
#     2022-Jul-28: move temperature control to function photosystem_temperature_dependence!
#     2023-Mar-11: add option to compute respiration rate only
#     2023-Sep-09: save ϕ_p to Leaf
#     2023-Oct-24: save ϕ_f1 and ϕ_f2
#
#######################################################################################################################################################################################################
# This method takes out stomtal model out and use it to determine whether to apply beta to Vcmax, Jmax, and Rd
leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT};
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(cache, leaf, air, leaf.flux.trait.stomatal_model; rd_only = rd_only);

# if stomtal model is not empirical model, then use the default β = 1
leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT},
            sm::AbstractStomataModel{FT};
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(cache, leaf, air, FT(1); rd_only = rd_only);

# if stomtal model is empirical model, then determine the β based on the parameter Y (if Vcmax, scale Vcmax, Jmax, and Rd)
leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT},
            sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}};
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(cache, leaf, air, sm.β, sm.β.PARAM_Y; rd_only = rd_only);

leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT},
            β::BetaFunction{FT},
            param_y::BetaParameterG1;
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(cache, leaf, air, FT(1); rd_only = rd_only);

leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT},
            β::BetaFunction{FT},
            param_y::BetaParameterVcmax;
            rd_only::Bool = false) where {FT} = leaf_photosynthesis!(cache, leaf, air, leaf.flux.auxil.β; rd_only = rd_only);

# This method computes and save the photosynthetic rates into leaf flux struct for Conductance mode
leaf_photosynthesis!(
            cache::SPACCache{FT},
            leaf::Leaf{FT},
            air::AirLayer{FT},
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

    # leaf._p_CO₂_i is not accurate here in the first call, thus need a second call after p_CO₂_i is analytically resolved
    # loop through sunlit leaf
    for i in eachindex(leaf.flux.auxil.ppar_sunlit)
        photosystem_electron_transport!(leaf.photosystem, leaf.flux.auxil.ppar_sunlit[i], leaf.flux.auxil.p_CO₂_i_sunlit[i]; β = β);
        rubisco_limited_rate!(cache, leaf.photosystem, air, leaf.flux.auxil.g_CO₂_sunlit[i]; β = β);
        light_limited_rate!(cache, leaf.photosystem, air, leaf.flux.auxil.g_CO₂_sunlit[i]; β = β);
        product_limited_rate!(leaf.photosystem, air, leaf.flux.auxil.g_CO₂_sunlit[i]; β = β);
        colimit_photosynthesis!(leaf.photosystem; β = β);

        # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
        leaf.flux.auxil.p_CO₂_i_sunlit[i] = air.s_aux.ps[2] - leaf.photosystem.auxil.a_n / leaf.flux.auxil.g_CO₂_sunlit[i] * air.state.p_air * FT(1e-6);
        leaf.flux.auxil.p_CO₂_s_sunlit[i] = air.s_aux.ps[2] - leaf.photosystem.auxil.a_n / leaf.flux.auxil.g_CO₂_b         * air.state.p_air * FT(1e-6);

        # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
        photosystem_electron_transport!(leaf.photosystem, leaf.flux.auxil.ppar_sunlit[i], leaf.flux.auxil.p_CO₂_i_sunlit[i]; β = β);

        # update the fluorescence related parameters
        photosystem_coefficients!(leaf.photosystem, leaf.flux.auxil.ppar_sunlit[i]; β = β);

        # save the rates and to leaf
        leaf.flux.auxil.a_n_sunlit[i]  = leaf.photosystem.auxil.a_n;
        leaf.flux.auxil.a_g_sunlit[i]  = leaf.photosystem.auxil.a_g;
        leaf.flux.auxil.etr_sunlit[i]  = leaf.photosystem.auxil.a_g / leaf.photosystem.auxil.e2c;
        leaf.flux.auxil.ϕ_f_sunlit[i]  = leaf.photosystem.auxil.ϕ_f;
        leaf.flux.auxil.ϕ_p_sunlit[i]  = leaf.photosystem.auxil.ϕ_p;
        leaf.flux.auxil.ϕ_f1_sunlit[i] = leaf.photosystem.auxil.ϕ_f1;
        leaf.flux.auxil.ϕ_f2_sunlit[i] = leaf.photosystem.auxil.ϕ_f2;
    end;

    # run the model for shaded leaf
    photosystem_electron_transport!(leaf.photosystem, leaf.flux.auxil.ppar_shaded, leaf.flux.auxil.p_CO₂_i_shaded; β = β);
    rubisco_limited_rate!(cache, leaf.photosystem, air, leaf.flux.auxil.g_CO₂_shaded; β = β);
    light_limited_rate!(cache, leaf.photosystem, air, leaf.flux.auxil.g_CO₂_shaded; β = β);
    product_limited_rate!(leaf.photosystem, air, leaf.flux.auxil.g_CO₂_shaded; β = β);
    colimit_photosynthesis!(leaf.photosystem; β = β);

    # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
    leaf.flux.auxil.p_CO₂_i_shaded = air.s_aux.ps[2] - leaf.photosystem.auxil.a_n / leaf.flux.auxil.g_CO₂_shaded * air.state.p_air * FT(1e-6);
    leaf.flux.auxil.p_CO₂_s_shaded = air.s_aux.ps[2] - leaf.photosystem.auxil.a_n / leaf.flux.auxil.g_CO₂_b      * air.state.p_air * FT(1e-6);

    # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
    photosystem_electron_transport!(leaf.photosystem, leaf.flux.auxil.ppar_shaded, leaf.flux.auxil.p_CO₂_i_shaded; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(leaf.photosystem, leaf.flux.auxil.ppar_shaded; β = β);

    # save the rates and to leaf
    leaf.flux.auxil.a_n_shaded = leaf.photosystem.auxil.a_n;
    leaf.flux.auxil.a_g_shaded = leaf.photosystem.auxil.a_g;
    leaf.flux.auxil.etr_shaded = leaf.photosystem.auxil.a_g / leaf.photosystem.auxil.e2c;
    leaf.flux.auxil.ϕ_f_shaded = leaf.photosystem.auxil.ϕ_f;
    leaf.flux.auxil.ϕ_p_shaded = leaf.photosystem.auxil.ϕ_p;

    return nothing
);
