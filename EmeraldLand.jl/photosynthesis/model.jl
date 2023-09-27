#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function leaf_photosynthesis!
#
#######################################################################################################################################################################################################
"""
Per refactored Photosynthesis module, the only things one need to know is the public function `leaf_photosynthesis!` and some construtors from `EmeraldNamespace`. See the examples in the methods
    below for details about how to use the function. The steps for computing photosynthetic rates are

- Update temperature dependent variables using [`photosystem_temperature_dependence!`](@ref)
- Calculate electron transport rate using [`photosystem_electron_transport!`](@ref)
- Calculate RubisCO limited rate using [`rubisco_limited_rate!`](@ref)
- Calculate light limited rate using [`light_limited_rate!`](@ref)
- Calculate product limited rate using [`product_limited_rate!`](@ref)
- Calculate gross and net rates using [`colimit_photosynthesis!`](@ref)
- Update fluorescence related variables using [`photosystem_coefficients!`](@ref)

"""
function leaf_photosynthesis! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method to compute photosynthetic rates only
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(lf::Union{Leaf2{FT}, Leaves2D{FT}}, air::AirLayer{FT}, g_lc::FT, ppar::FT, t::FT = lf.t) where {FT}

Updates leaf photosynthetic rates based on CO₂ partial pressure (for StomataModels.jl temporary use), given
- `lf` `Leaf2`, `Leaves2D` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, default is `leaf._g_CO₂`
- `ppar` APAR used for photosynthesis
- `t` Leaf temperature in `[K]`

"""
leaf_photosynthesis!(lf::Union{Leaf2{FT}, Leaves2D{FT}}, air::AirLayer{FT}, g_lc::FT, ppar::FT, t::FT = lf.t) where {FT} = (
    (; PRC, PSM) = lf;

    photosystem_temperature_dependence!(PSM, PRC, air, t);
    photosystem_electron_transport!(PSM, PRC, ppar, FT(20); β = FT(1));
    rubisco_limited_rate!(PSM, air, g_lc; β = FT(1));
    light_limited_rate!(PSM, PRC, air, g_lc; β = FT(1));
    product_limited_rate!(PSM, air, g_lc; β = FT(1));
    colimit_photosynthesis!(PSM; β = FT(1));

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-12: add method to account for tuning factor at leaf level
#     2023-Mar-11: add option to compute respiration rate only
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(lf::Union{Leaf{FT}, Leaves2D{FT}}, air::AirLayer{FT}, mode::Union{GCO₂Mode, PCO₂Mode}; rd_only::Bool = false) where {FT}

Updates leaf photosynthetic rates based on CO₂ partial pressure or CO₂ conductance, given
- `lf` `Leaf`, `Leaves2D` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `GCO₂Mode` or `PCO₂Mode` that uses CO₂ conductance or partial pressure to compute photosynthetic rates
- `rd_only` Whether to compute respiration rate only

"""
leaf_photosynthesis!(lf::Union{Leaf2{FT}, Leaves2D{FT}}, air::AirLayer{FT}, mode::Union{GCO₂Mode, PCO₂Mode}; rd_only::Bool = false) where {FT} =
    leaf_photosynthesis!(lf, air, mode, lf.SM; rd_only = rd_only);

leaf_photosynthesis!(
            lf::Union{Leaf2{FT}, Leaves2D{FT}},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            sm::AbstractStomataModel{FT};
            rd_only::Bool = false
) where {FT} = leaf_photosynthesis!(lf, air, mode, FT(1); rd_only = rd_only);

leaf_photosynthesis!(
            lf::Union{Leaf2{FT}, Leaves2D{FT}},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}};
            rd_only::Bool = false
) where {FT} = leaf_photosynthesis!(lf, air, mode, sm.β, sm.β.PARAM_Y; rd_only = rd_only);

leaf_photosynthesis!(
            lf::Union{Leaf2{FT}, Leaves2D{FT}},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            β::BetaFunction{FT},
            param_y::BetaParameterG1;
            rd_only::Bool = false
) where {FT} = leaf_photosynthesis!(lf, air, mode, FT(1); rd_only = rd_only);

leaf_photosynthesis!(
            lf::Union{Leaf2{FT}, Leaves2D{FT}},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            β::BetaFunction{FT},
            param_y::BetaParameterVcmax;
            rd_only::Bool = false
) where {FT} = leaf_photosynthesis!(lf, air, mode, β.β₁; rd_only = rd_only);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: set a default p_i from leaf to combine two methods
#     2022-Jan-14: do not update temperature to avoid its impact on plant hydraulics
#     2022-Jan-14: use colimit function to compute a_gross and a_net
#     2022-Jan-14: set a default g_lc from leaf to combine two methods
#     2022-Jan-18: add p_i to electron transport function input variables
#     2022-Jan-24: fix PSM abstraction in colimit_photosynthesis! function
#     2022-Feb-07: use new method of photosystem_coefficients!
#     2022-Feb-28: use updated light_limited_rate! function
#     2022-Feb-28: add support to C3CytochromeModel
#     2022-Jun-27: remove ppar from input variable list of light_limited_rate!
#     2022-Jun-28: add method for Leaves2D
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2022-Jul-07: save a_net and a_gross to Leaf (as PSM may be used for temporary calculations)
#     2022-Jul-12: use β as a must have option (and thus this function becomes a core function of the one above)
#     2022-Jul-28: move temperature control to function photosystem_temperature_dependence!
#     2023-Mar-11: add option to compute respiration rate only
#     2023-Jun-13: save actual etr as well
#     2023-Sep-09: save ϕ_d, ϕ_n, and ϕ_p to Leaves2D
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(leaf::Leaf2{FT}, air::AirLayer{FT}, mode::PCO₂Mode, β::FT; rd_only::Bool = false) where {FT}
    leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::PCO₂Mode, β::FT; rd_only::Bool = false) where {FT}
    leaf_photosynthesis!(leaf::Leaf2{FT}, air::AirLayer{FT}, mode::GCO₂Mode, β::FT; rd_only::Bool = false) where {FT}
    leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::GCO₂Mode, β::FT; rd_only::Bool = false) where {FT}

Updates leaf photosynthetic rates (this method not meant for public usage, use it with caution), given
- `leaf` `Leaf2` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `leaves` `Leaves2D` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `GCO₂Mode` or `PCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd
- `rd_only` Whether to compute respiration rate only

"""
leaf_photosynthesis!(leaf::Leaf2{FT}, air::AirLayer{FT}, mode::PCO₂Mode, β::FT; rd_only::Bool = false) where {FT} = (
    (; PRC, PSM) = leaf;

    if rd_only
        PSM._r_d = PSM.r_d25 * temperature_correction(PSM.TD_R, leaf.t);
        leaf.a_net = -PSM._r_d;
        leaf.a_gross = 0;
        leaf.etr = 0;

        return nothing
    end;

    photosystem_temperature_dependence!(PSM, PRC, air, leaf.t);
    photosystem_electron_transport!(PSM, PRC, leaf.ppar, leaf._p_CO₂_i; β = β);
    rubisco_limited_rate!(PSM, leaf._p_CO₂_i; β = β);
    light_limited_rate!(PSM);
    product_limited_rate!(PSM, leaf._p_CO₂_i; β = β);
    colimit_photosynthesis!(PSM; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaf.ppar; β = β);

    # save the rates and to leaf
    leaf.a_net = PSM.a_net;
    leaf.a_gross = PSM.a_gross;
    leaf.etr = PSM.a_gross / PSM._e_to_c;

    return nothing
);

leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::PCO₂Mode, β::FT; rd_only::Bool = false) where {FT} = (
    (; PRC, PSM) = leaves;

    if rd_only
        PSM._r_d = PSM.r_d25 * temperature_correction(PSM.TD_R, leaves.t);
        leaves.a_net_sunlit .= -PSM._r_d;
        leaves.a_gross_sunlit .= 0;
        leaves.etr_sunlit .= 0;
        leaves.ϕ_f_sunlit .= 0;
        leaves.a_net_shaded = -PSM._r_d;
        leaves.a_gross_shaded = 0;
        leaves.etr_shaded = 0;
        leaves.ϕ_f_shaded = 0;

        return nothing
    end;

    photosystem_temperature_dependence!(PSM, PRC, air, leaves.t);

    # loop through the ppars for sunlit leaves
    for _i in eachindex(leaves.ppar_sunlit)
        # calculate the photosynthetic rates
        photosystem_electron_transport!(PSM, PRC, leaves.ppar_sunlit[_i], leaves._p_CO₂_i_sunlit[_i]; β = β);
        rubisco_limited_rate!(PSM, leaves._p_CO₂_i_sunlit[_i]; β = β);
        light_limited_rate!(PSM);
        product_limited_rate!(PSM, leaves._p_CO₂_i_sunlit[_i]; β = β);
        colimit_photosynthesis!(PSM; β = β);

        # update the fluorescence related parameters
        photosystem_coefficients!(PSM, PRC, leaves.ppar_sunlit[_i]; β = β);

        # save the rates and to leaves
        leaves.a_net_sunlit[_i] = PSM.a_net;
        leaves.a_gross_sunlit[_i] = PSM.a_gross;
        leaves.etr_sunlit[_i] = PSM.a_gross / PSM._e_to_c;
        leaves.ϕ_d_sunlit[_i] = PRC.ϕ_d;
        leaves.ϕ_f_sunlit[_i] = PRC.ϕ_f;
        leaves.ϕ_n_sunlit[_i] = PRC.ϕ_n;
        leaves.ϕ_p_sunlit[_i] = PRC.ϕ_p;
    end;

    # run the model for shaded leaves
    photosystem_electron_transport!(PSM, PRC, leaves.ppar_shaded, leaves._p_CO₂_i_shaded; β = β);
    rubisco_limited_rate!(PSM, leaves._p_CO₂_i_shaded; β = β);
    light_limited_rate!(PSM);
    product_limited_rate!(PSM, leaves._p_CO₂_i_shaded; β = β);
    colimit_photosynthesis!(PSM; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaves.ppar_shaded; β = β);

    # save the rates and to leaves
    leaves.a_net_shaded = PSM.a_net;
    leaves.a_gross_shaded = PSM.a_gross;
    leaves.etr_shaded = PSM.a_gross / PSM._e_to_c;
    leaves.ϕ_d_shaded = PRC.ϕ_d;
    leaves.ϕ_f_shaded = PRC.ϕ_f;
    leaves.ϕ_n_shaded = PRC.ϕ_n;
    leaves.ϕ_p_shaded = PRC.ϕ_p;

    return nothing
);

leaf_photosynthesis!(leaf::Leaf2{FT}, air::AirLayer{FT}, mode::GCO₂Mode, β::FT; rd_only::Bool = false) where {FT} = (
    (; PRC, PSM) = leaf;

    if rd_only
        PSM._r_d = PSM.r_d25 * temperature_correction(PSM.TD_R, leaf.t);
        leaf.a_net = -PSM._r_d;
        leaf.a_gross = 0;
        leaf.etr = 0;

        return nothing
    end;

    # leaf._p_CO₂_i is not accurate here in the first call, thus need a second call after p_CO₂_i is analytically resolved
    photosystem_temperature_dependence!(PSM, PRC, air, leaf.t);
    photosystem_electron_transport!(PSM, PRC, leaf.ppar, leaf._p_CO₂_i; β = β);
    rubisco_limited_rate!(PSM, air, leaf._g_CO₂; β = β);
    light_limited_rate!(PSM, PRC, air, leaf._g_CO₂; β = β);
    product_limited_rate!(PSM, air, leaf._g_CO₂; β = β);
    colimit_photosynthesis!(PSM; β = β);

    # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
    leaf._p_CO₂_i = air.p_CO₂ - PSM.a_net / leaf._g_CO₂   * air.P_AIR * FT(1e-6);
    leaf._p_CO₂_s = air.p_CO₂ - PSM.a_net / leaf.g_CO₂_b * air.P_AIR * FT(1e-6);

    # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
    photosystem_electron_transport!(PSM, PRC, leaf.ppar, leaf._p_CO₂_i; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaf.ppar; β = β);

    # save the rates and to leaf
    leaf.a_net = PSM.a_net;
    leaf.a_gross = PSM.a_gross;
    leaf.etr = PSM.a_gross / PSM._e_to_c;

    return nothing
);

leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::GCO₂Mode, β::FT; rd_only::Bool = false) where {FT} = (
    (; PRC, PSM) = leaves;

    if rd_only
        PSM._r_d = PSM.r_d25 * temperature_correction(PSM.TD_R, leaves.t);
        leaves.a_net_sunlit .= -PSM._r_d;
        leaves.a_gross_sunlit .= 0;
        leaves.etr_sunlit .= 0;
        leaves.ϕ_f_sunlit .= 0;
        leaves.a_net_shaded = -PSM._r_d;
        leaves.a_gross_shaded = 0;
        leaves.etr_shaded = 0;
        leaves.ϕ_f_shaded = 0;

        return nothing
    end;

    photosystem_temperature_dependence!(PSM, PRC, air, leaves.t);

    # leaf._p_CO₂_i is not accurate here in the first call, thus need a second call after p_CO₂_i is analytically resolved
    # loop through sunlit leaves
    for _i in eachindex(leaves.ppar_sunlit)
        photosystem_electron_transport!(PSM, PRC, leaves.ppar_sunlit[_i], leaves._p_CO₂_i_sunlit[_i]; β = β);
        rubisco_limited_rate!(PSM, air, leaves._g_CO₂_sunlit[_i]; β = β);
        light_limited_rate!(PSM, PRC, air, leaves._g_CO₂_sunlit[_i]; β = β);
        product_limited_rate!(PSM, air, leaves._g_CO₂_sunlit[_i]; β = β);
        colimit_photosynthesis!(PSM; β = β);

        # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
        leaves._p_CO₂_i_sunlit[_i] = air.p_CO₂ - PSM.a_net / leaves._g_CO₂_sunlit[_i] * air.P_AIR * FT(1e-6);
        leaves._p_CO₂_s_sunlit[_i] = air.p_CO₂ - PSM.a_net / leaves.g_CO₂_b          * air.P_AIR * FT(1e-6);

        # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
        photosystem_electron_transport!(PSM, PRC, leaves.ppar_sunlit[_i], leaves._p_CO₂_i_sunlit[_i]; β = β);

        # update the fluorescence related parameters
        photosystem_coefficients!(PSM, PRC, leaves.ppar_sunlit[_i]; β = β);

        # save the rates and to leaves
        leaves.a_net_sunlit[_i] = PSM.a_net;
        leaves.a_gross_sunlit[_i] = PSM.a_gross;
        leaves.etr_sunlit[_i] = PSM.a_gross / PSM._e_to_c;
        leaves.ϕ_d_sunlit[_i] = PRC.ϕ_d;
        leaves.ϕ_f_sunlit[_i] = PRC.ϕ_f;
        leaves.ϕ_n_sunlit[_i] = PRC.ϕ_n;
        leaves.ϕ_p_sunlit[_i] = PRC.ϕ_p;
    end;

    # run the model for shaded leaves
    photosystem_electron_transport!(PSM, PRC, leaves.ppar_shaded, leaves._p_CO₂_i_shaded; β = β);
    rubisco_limited_rate!(PSM, air, leaves._g_CO₂_shaded; β = β);
    light_limited_rate!(PSM, PRC, air, leaves._g_CO₂_shaded; β = β);
    product_limited_rate!(PSM, air, leaves._g_CO₂_shaded; β = β);
    colimit_photosynthesis!(PSM; β = β);

    # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
    leaves._p_CO₂_i_shaded = air.p_CO₂ - PSM.a_net / leaves._g_CO₂_shaded * air.P_AIR * FT(1e-6);
    leaves._p_CO₂_s_shaded = air.p_CO₂ - PSM.a_net / leaves.g_CO₂_b      * air.P_AIR * FT(1e-6);

    # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
    photosystem_electron_transport!(PSM, PRC, leaves.ppar_shaded, leaves._p_CO₂_i_shaded; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaves.ppar_shaded; β = β);

    # save the rates and to leaves
    leaves.a_net_shaded = PSM.a_net;
    leaves.a_gross_shaded = PSM.a_gross;
    leaves.etr_shaded = PSM.a_gross / PSM._e_to_c;
    leaves.ϕ_d_shaded = PRC.ϕ_d;
    leaves.ϕ_f_shaded = PRC.ϕ_f;
    leaves.ϕ_n_shaded = PRC.ϕ_n;
    leaves.ϕ_p_shaded = PRC.ϕ_p;

    return nothing
);


######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-29: add method for MonoElementSPAC
#     2022-Jun-29: add method for MultiLayerSPAC
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2022-Jul-13: redirect the wrapper function to the method at leaf level
#     2023-Mar-11: only compute respiration rate if solar zenith angle >= 89
#     2023-Mar-11: do nothing if LAI == 0
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(spac::MonoElementSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT}
    leaf_photosynthesis!(spac::MultiLayerSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT}

Updates leaf photosynthetic rates for SPAC, given
- `spac` `MonoElementSPAC` or `MultiLayerSPAC` type SPAC
- `mode` `GCO₂Mode` or `PCO₂Mode`

"""
leaf_photosynthesis!(spac::MonoElementSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT} = leaf_photosynthesis!(spac.LEAF, spac.AIR, mode);

leaf_photosynthesis!(spac::MultiLayerSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT} = (
    (; AIR, ANGLES, CANOPY, LEAVES, LEAVES_INDEX) = spac;

    if CANOPY.lai == 0
        return nothing
    end;

    _rd_only = ANGLES.sza < 89 ? false : true;
    for _i in eachindex(LEAVES)
        leaf_photosynthesis!(LEAVES[_i], AIR[LEAVES_INDEX[_i]], mode; rd_only = _rd_only);
    end;

    return nothing
);
