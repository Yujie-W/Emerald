#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Apr-13: name the method to shortwave_radiation! to be more specific
#
#######################################################################################################################################################################################################
"""

Run short wave radiation.

"""
function shortwave_radiation! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-14: make method work with broadband soil albedo struct
#     2022-Jun-14: allow method to use broadband PAR and NIR soil albedo values
#     2023-Apr-13: name the method to shortwave_radiation! to be more specific
#
#######################################################################################################################################################################################################
"""

    shortwave_radiation!(can::HyperspectralMLCanopy{FT}, albedo::BroadbandSoilAlbedo{FT}) where {FT}
    shortwave_radiation!(can::HyperspectralMLCanopy{FT}, albedo::HyperspectralSoilAlbedo{FT}) where {FT}

Updates soil shortwave radiation profiles, given
- `config` Configuration for `MultiLayerSPAC`
- `can` `HyperspectralMLCanopy` type struct
- `albedo` `BroadbandSoilAlbedo` or `HyperspectralSoilAlbedo` type soil albedo

"""
shortwave_radiation!(config::SPACConfiguration{FT}, can::HyperspectralMLCanopy{FT}, albedo::BroadbandSoilAlbedo{FT}) where {FT} = (
    (; SPECTRA) = config;
    (; DIM_LAYER, OPTICS, RADIATION) = can;

    OPTICS._tmp_vec_λ[SPECTRA.IΛ_PAR] .= view(RADIATION.e_direct,SPECTRA.IΛ_PAR,DIM_LAYER+1) .* (1 .- albedo.ρ_sw[1]);
    OPTICS._tmp_vec_λ[SPECTRA.IΛ_NIR] .= view(RADIATION.e_direct,SPECTRA.IΛ_NIR,DIM_LAYER+1) .* (1 .- albedo.ρ_sw[2]);
    albedo.e_net_direct = OPTICS._tmp_vec_λ' * SPECTRA.ΔΛ / 1000;

    OPTICS._tmp_vec_λ[SPECTRA.IΛ_PAR] .= view(RADIATION.e_diffuse_down,SPECTRA.IΛ_PAR,DIM_LAYER+1) .* (1 .- albedo.ρ_sw[1]);
    OPTICS._tmp_vec_λ[SPECTRA.IΛ_NIR] .= view(RADIATION.e_diffuse_down,SPECTRA.IΛ_NIR,DIM_LAYER+1) .* (1 .- albedo.ρ_sw[2]);
    albedo.e_net_diffuse = OPTICS._tmp_vec_λ' * SPECTRA.ΔΛ / 1000;

    albedo.r_net_sw = albedo.e_net_direct + albedo.e_net_diffuse;

    return nothing
);

shortwave_radiation!(config::SPACConfiguration{FT}, can::HyperspectralMLCanopy{FT}, albedo::HyperspectralSoilAlbedo{FT}) where {FT} = (
    (; DIM_LAYER, SPECTRA) = config;
    (; RADIATION) = can;

    albedo.e_net_direct .= view(RADIATION.e_direct,:,DIM_LAYER+1) .* (1 .- albedo.ρ_sw);
    albedo.e_net_diffuse .= view(RADIATION.e_diffuse_down,:,DIM_LAYER+1) .* (1 .- albedo.ρ_sw);
    albedo.r_net_sw = (albedo.e_net_direct' * SPECTRA.ΔΛ + albedo.e_net_diffuse' * SPECTRA.ΔΛ) / 1000;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-09: migrate the function from CanopyLayers
#     2022-Jun-10: rename PAR/APAR to APAR/PPAR to be more accurate
#     2022-Jun-10: add PAR calculation (before absorption)
#     2022-Jun-10: compute shortwave net radiation
#     2022-Jun-13: use DIM_LAYER instead of _end
#     2022-Jun-29: use Leaf for the hyperspectral RT
#     2023-Mar-11: add code to account for the case of LAI == 0
#     2023-Apr-13: rename option APAR_car to apar_car
#     2023-Apr-13: name the method to shortwave_radiation! to be more specific
#     2023-May-19: use δlai per canopy layer
#     2023-Jun-15: compute PAR when lai = 0
#     2023-Jun-20: remove option apar_car as it is already in config
# Bug fixes
#     2022-Jul-15: sum by r_net_sw by the weights of sunlit and shaded fractions
#     2022-Jul-27: use _ρ_dd, _ρ_sd, _τ_dd, and _τ_sd for leaf energy absorption (typo when refactoring the code)
#     2022-Aug-30: fix par, apar, and ppar issues
#     2022-Nov-29: ppar for sunlit leaves should be scaled based on sunlit fraction
#
#######################################################################################################################################################################################################
"""

    shortwave_radiation!(config::SPACConfiguration{FT}, can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::HyperspectralRadiation{FT}, soil::Soil{FT}) where {FT}

Updates canopy radiation profiles for shortwave radiation, given
- `config` Configuration for `MultiLayerSPAC`
- `can` `HyperspectralMLCanopy` type struct
- `leaves` Vector of `Leaf`
- `rad` Incoming shortwave radiation
- `soil` Bottom soil boundary layer

"""
shortwave_radiation!(config::SPACConfiguration{FT}, can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::HyperspectralRadiation{FT}, soil::Soil{FT}) where {FT} = (
    (; DIM_LAYER, SPECTRA) = config;
    (; OPTICS, P_INCL, RADIATION) = can;
    (; ALBEDO) = soil;

    if can.lai == 0
        for _i in 1:(DIM_LAYER+1)
            RADIATION.e_direct[:,_i] .= rad.e_direct;
            RADIATION.e_diffuse_down[:,_i] .= rad.e_diffuse;
            RADIATION.e_diffuse_up[:,_i] .= 0;
            RADIATION.e_v[:,_i] .= 0;
        end;
        RADIATION.e_diffuse_up[:,end] .= view(OPTICS.ρ_sd,:,DIM_LAYER+1) .* view(RADIATION.e_direct,:,DIM_LAYER+1) .+ view(OPTICS.ρ_dd,:,DIM_LAYER+1) .* view(RADIATION.e_diffuse_down,:,DIM_LAYER+1);
        RADIATION.e_v[:,end] .= view(RADIATION.e_diffuse_up,:,DIM_LAYER+1);
        RADIATION.e_o .= view(RADIATION.e_diffuse_up,:,DIM_LAYER+1) ./ FT(pi);
        RADIATION.albedo .= RADIATION.e_o * FT(pi) ./ (rad.e_direct .+ rad.e_diffuse);

        RADIATION._par_shaded .= photon.(SPECTRA.Λ_PAR, view(rad.e_diffuse,SPECTRA.IΛ_PAR)) .* 1000;
        RADIATION._par_sunlit .= photon.(SPECTRA.Λ_PAR, view(rad.e_direct ,SPECTRA.IΛ_PAR)) .* 1000;
        RADIATION.par_in_diffuse = RADIATION._par_shaded' * SPECTRA.ΔΛ_PAR;
        RADIATION.par_in_direct = RADIATION._par_sunlit' * SPECTRA.ΔΛ_PAR;
        RADIATION.par_in = RADIATION.par_in_diffuse + RADIATION.par_in_direct;

        return nothing
    end;

    # 1. update upward and downward direct and diffuse radiation profiles
    RADIATION.e_direct[:,1] .= rad.e_direct;
    RADIATION.e_diffuse_down[:,1] .= rad.e_diffuse;

    for _i in 1:DIM_LAYER
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i  );     # downward diffuse radiation at upper boundary
        _e_d_j = view(RADIATION.e_diffuse_down,:,_i+1);     # downward diffuse radiation at lower boundary
        _e_s_i = view(RADIATION.e_direct      ,:,_i  );     # direct radiation at upper boundary
        _e_s_j = view(RADIATION.e_direct      ,:,_i+1);     # direct radiation at lower boundary
        _e_u_i = view(RADIATION.e_diffuse_up  ,:,_i  );     # upward diffuse radiation at upper boundary

        _r_dd_i = view(OPTICS.ρ_dd ,:,_i);  # reflectance of the upper boundary (i)
        _r_sd_i = view(OPTICS.ρ_sd ,:,_i);  # reflectance of the upper boundary (i)
        _t_dd_i = view(OPTICS.τ_dd ,:,_i);  # transmittance of the layer (i)
        _t_sd_i = view(OPTICS.τ_sd ,:,_i);  # transmittance of the layer (i)
        _t_ss__ = view(OPTICS._τ_ss,_i);    # transmittance for directional->directional

        _e_s_j .= _t_ss__ .* _e_s_i;
        _e_d_j .= _t_sd_i .* _e_s_i .+ _t_dd_i .* _e_d_i;
        _e_u_i .= _r_sd_i .* _e_s_i .+ _r_dd_i .* _e_d_i;
    end;

    RADIATION.e_diffuse_up[:,end] .= view(OPTICS.ρ_sd,:,DIM_LAYER+1) .* view(RADIATION.e_direct,:,DIM_LAYER+1) .+ view(OPTICS.ρ_dd,:,DIM_LAYER+1) .* view(RADIATION.e_diffuse_down,:,DIM_LAYER+1);

    # 2. update the sunlit and shaded sum radiation and total absorbed radiation per layer and for soil
    for _i in 1:DIM_LAYER
        _a_s_i = view(RADIATION.e_net_direct  ,:,_i  );     # net absorbed direct radiation
        _a_d_i = view(RADIATION.e_net_diffuse ,:,_i  );     # net absorbed diffuse radiation
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i  );     # downward diffuse radiation at upper boundary
        _e_s_i = view(RADIATION.e_direct      ,:,_i  );     # direct radiation at upper boundary
        _e_u_j = view(RADIATION.e_diffuse_up  ,:,_i+1);     # upward diffuse radiation at lower boundary
        _p_s_i = view(RADIATION.e_sum_direct  ,:,_i  );     # sum direct radiation
        _p_d_i = view(RADIATION.e_sum_diffuse ,:,_i  );     # sum diffuse radiation

        _r_dd__ = view(OPTICS._ρ_dd,:,_i);  # reflectance of the upper boundary (i)
        _r_sd__ = view(OPTICS._ρ_sd,:,_i);  # reflectance of the upper boundary (i)
        _t_dd__ = view(OPTICS._τ_dd,:,_i);  # transmittance of the layer (i)
        _t_sd__ = view(OPTICS._τ_sd,:,_i);  # transmittance of the layer (i)
        _t_ss__ = view(OPTICS._τ_ss,_i);    # transmittance for directional->directional

        _p_s_i .= _e_s_i;
        _p_d_i .= _e_d_i .+ _e_u_j;

        _a_s_i .= _p_s_i .* (1 .- _t_ss__ .- _t_sd__ .- _r_sd__);
        _a_d_i .= _p_d_i .* (1 .- _t_dd__ .- _r_dd__);
    end;

    # 3. compute the spectra at the observer direction
    for _i in 1:DIM_LAYER
        _e_d_i = view(RADIATION.e_diffuse_down,:,_i);   # downward diffuse radiation at upper boundary
        _e_u_i = view(RADIATION.e_diffuse_up  ,:,_i);   # upward diffuse radiation at upper boundary
        _e_v_i = view(RADIATION.e_v,:,_i);

        _dob_i = view(OPTICS.σ_dob,:,_i);   # scattering coefficient backward for diffuse->observer
        _dof_i = view(OPTICS.σ_dof,:,_i);   # scattering coefficient forward for diffuse->observer
        _so__i = view(OPTICS.σ_so ,:,_i);   # bidirectional from solar to observer

        # _e_v_i .= (OPTICS.po[_i] .* _dob_i .* _e_d_i .+ OPTICS.po[_i] .* _dof_i .* _e_u_i .+ OPTICS.pso[_i] .* _so__i .* rad.e_direct) * _ilai;
        _e_v_i  .= OPTICS.po[_i] .* _dob_i .* _e_d_i;
        _e_v_i .+= OPTICS.po[_i] .* _dof_i .* _e_u_i;
        _e_v_i .+= OPTICS.pso[_i] .* _so__i .* rad.e_direct;
        _e_v_i .*= can.δlai[_i] * can.ci;
    end;
    RADIATION.e_v[:,end] .= OPTICS.po[end] .* view(RADIATION.e_diffuse_up,:,DIM_LAYER+1);

    for _i in eachindex(RADIATION.e_o)
        RADIATION.e_o[_i] = sum(view(RADIATION.e_v,_i,:)) / FT(pi);
    end;

    RADIATION.albedo .= RADIATION.e_o * FT(pi) ./ (rad.e_direct .+ rad.e_diffuse);

    # 4. compute net absorption for leaves and soil
    for _i in 1:DIM_LAYER
        _Σ_shaded = view(RADIATION.e_net_diffuse,:,_i)' * SPECTRA.ΔΛ / 1000 / can.δlai[_i];
        _Σ_sunlit = view(RADIATION.e_net_direct ,:,_i)' * SPECTRA.ΔΛ / 1000 / can.δlai[_i];
        RADIATION.r_net_sw_shaded[_i] = _Σ_shaded;
        RADIATION.r_net_sw_sunlit[_i] = _Σ_sunlit / OPTICS.p_sunlit[_i] + _Σ_shaded;
        RADIATION.r_net_sw[_i] = _Σ_shaded * (1 - OPTICS.p_sunlit[_i]) + _Σ_sunlit * OPTICS.p_sunlit[_i];
    end;

    shortwave_radiation!(config, can, ALBEDO);

    # 5. compute top-of-canopy and leaf level PAR, APAR, and PPAR per ground area
    RADIATION._par_shaded .= photon.(SPECTRA.Λ_PAR, view(rad.e_diffuse,SPECTRA.IΛ_PAR)) .* 1000;
    RADIATION._par_sunlit .= photon.(SPECTRA.Λ_PAR, view(rad.e_direct ,SPECTRA.IΛ_PAR)) .* 1000;
    RADIATION.par_in_diffuse = RADIATION._par_shaded' * SPECTRA.ΔΛ_PAR;
    RADIATION.par_in_direct = RADIATION._par_sunlit' * SPECTRA.ΔΛ_PAR;
    RADIATION.par_in = RADIATION.par_in_diffuse + RADIATION.par_in_direct;

    mul!(OPTICS._tmp_vec_azi, OPTICS._abs_fs', P_INCL);
    _normi = 1 / mean(OPTICS._tmp_vec_azi);

    for _i in 1:DIM_LAYER
        _α_apar = view(leaves[_i].bio.auxil.f_ppar, SPECTRA.IΛ_PAR);

        # convert energy to quantum unit for PAR, APAR and PPAR per leaf area
        RADIATION._par_shaded  .= photon.(SPECTRA.Λ_PAR, view(RADIATION.e_sum_diffuse,SPECTRA.IΛ_PAR,_i)) .* 1000;
        RADIATION._par_sunlit  .= photon.(SPECTRA.Λ_PAR, view(RADIATION.e_sum_direct ,SPECTRA.IΛ_PAR,_i)) .* 1000 ./ OPTICS.p_sunlit[_i];
        RADIATION._apar_shaded .= photon.(SPECTRA.Λ_PAR, view(RADIATION.e_net_diffuse,SPECTRA.IΛ_PAR,_i)) .* 1000 ./ can.δlai[_i];
        RADIATION._apar_sunlit .= photon.(SPECTRA.Λ_PAR, view(RADIATION.e_net_direct ,SPECTRA.IΛ_PAR,_i)) .* 1000 ./ can.δlai[_i] ./ OPTICS.p_sunlit[_i];
        RADIATION._ppar_shaded .= RADIATION._apar_shaded .* _α_apar;
        RADIATION._ppar_sunlit .= RADIATION._apar_sunlit .* _α_apar;

        # PAR for leaves
        _Σ_par_dif = RADIATION._par_shaded' * SPECTRA.ΔΛ_PAR;
        _Σ_par_dir = RADIATION._par_sunlit' * SPECTRA.ΔΛ_PAR * _normi;
        RADIATION.par_shaded[_i] = _Σ_par_dif;
        RADIATION.par_sunlit[:,:,_i] .= OPTICS._abs_fs_fo .* _Σ_par_dir;
        RADIATION.par_sunlit[:,:,_i] .+= _Σ_par_dif;

        # APAR for leaves
        _Σ_apar_dif = RADIATION._apar_shaded' * SPECTRA.ΔΛ_PAR;
        _Σ_apar_dir = RADIATION._apar_sunlit' * SPECTRA.ΔΛ_PAR * _normi;
        RADIATION.apar_shaded[_i] = _Σ_apar_dif;
        RADIATION.apar_sunlit[:,:,_i] .= OPTICS._abs_fs_fo .* _Σ_apar_dir;
        RADIATION.apar_sunlit[:,:,_i] .+= _Σ_apar_dif;

        # PPAR for leaves
        _Σ_ppar_dif = RADIATION._ppar_shaded' * SPECTRA.ΔΛ_PAR;
        _Σ_ppar_dir = RADIATION._ppar_sunlit' * SPECTRA.ΔΛ_PAR * _normi;
        leaves[DIM_LAYER+1-_i].ppar_shaded  = _Σ_ppar_dif;
        leaves[DIM_LAYER+1-_i].ppar_sunlit .= OPTICS._abs_fs_fo .* _Σ_ppar_dir .+ _Σ_ppar_dif;
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Apr-13: move method of canopy_radiation! to longwave_radiation!
#
#######################################################################################################################################################################################################
"""

Update long wave radiation profiles for canopy.

"""
function longwave_radiation! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-10: migrate the function thermal_fluxes! from CanopyLayers
#     2022-Jun-10: update net lw radiation for leaves and soil
#     2022-Jun-29: use Leaf for the hyperspectral RT
#     2023-Mar-11: add code to account for the case of LAI == 0
#
#######################################################################################################################################################################################################
"""

    longwave_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::FT, soil::Soil{FT}) where {FT}

Updates canopy radiation profiles for shortwave or longwave radiation, given
- `can` `HyperspectralMLCanopy` type struct
- `leaves` Vector of `Leaf`
- `rad` Incoming longwave radiation
- `soil` Bottom soil boundary layer

"""
longwave_radiation!(can::HyperspectralMLCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::FT, soil::Soil{FT}) where {FT} = (
    (; OPTICS, RADIATION) = can;
    (; ALBEDO, LAYERS) = soil;

    _nlayers = length(leaves);

    if can.lai == 0
        _r_lw_soil = K_STEFAN(FT) * (1 - ALBEDO.ρ_LW) * LAYERS[1].t ^ 4;
        RADIATION.r_lw .= 0;
        RADIATION.r_net_lw .= 0;
        RADIATION.r_lw_up .= rad * ALBEDO.ρ_LW + _r_lw_soil;
        ALBEDO.r_net_lw = rad * (1 - ALBEDO.ρ_LW) - _r_lw_soil;

        return nothing
    end;

    # 1. compute longwave radiation out from the leaves and soil
    for _i in eachindex(leaves)
        RADIATION.r_lw[_i] = K_STEFAN(FT) * OPTICS.ϵ[_i] * leaves[_i].energy.auxil.t ^ 4;
    end;

    _r_lw_soil = K_STEFAN(FT) * (1 - ALBEDO.ρ_LW) * LAYERS[1].t ^ 4;

    # 2. account for the longwave emission from bottom to up
    RADIATION._r_emit_up[end] = _r_lw_soil;

    for _i in _nlayers:-1:1
        _r__ = OPTICS._ρ_lw[_i];
        _r_j = OPTICS.ρ_lw[_i+1];
        _t__ = OPTICS._τ_lw[_i];

        _dnorm = 1 - _r__ * _r_j;

        RADIATION._r_emit_down[_i] = (RADIATION._r_emit_up[_i+1] * _r__ + RADIATION.r_lw[_i]) / _dnorm;
        RADIATION._r_emit_up[_i] = RADIATION._r_emit_down[_i] * _r_j * _t__ + RADIATION._r_emit_up[_i+1] * _t__ + RADIATION.r_lw[_i];
    end;

    # 3. account for the longwave emission from up to bottom
    RADIATION.r_lw_down[1] = rad;

    for _i in 1:_nlayers
        _r_i = OPTICS.ρ_lw[_i];
        _t_i = OPTICS.τ_lw[_i];

        RADIATION.r_lw_down[_i+1] = RADIATION.r_lw_down[_i] * _t_i + RADIATION._r_emit_down[_i];
        RADIATION.r_lw_up[_i+1] = RADIATION.r_lw_down[_i] * _r_i + RADIATION._r_emit_up[_i];
    end;

    RADIATION.r_lw_up[end] = RADIATION.r_lw_down[end] * ALBEDO.ρ_LW + _r_lw_soil;

    # 4. compute the net longwave radiation per canopy layer and soil
    for _i in 1:_nlayers
        RADIATION.r_net_lw[_i] = (RADIATION.r_lw_down[_i] + RADIATION.r_lw_up[_i+1]) * (1 - OPTICS._ρ_lw[_i] - OPTICS._τ_lw[_i]) - 2 * RADIATION.r_lw[_i];
    end;

    ALBEDO.r_net_lw = RADIATION.r_lw_down[end] * (1 - ALBEDO.ρ_LW) - _r_lw_soil;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-09: migrate the function from CanopyLayers
#     2022-Jun-09: rename function to canopy_radiation!
#
#######################################################################################################################################################################################################
"""

Run shortwave and longwave radiation.

"""
function canopy_radiation! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-29: add method for SPAC
#     2022-Jul-28: update soil albedo at the very first step
#     2023-Mar-11: run canopy optical properties and shortwave radiation only if solar zenith angle is lower than 89
#     2023-Apr-13: sw and lw radiation moved to METEO
#     2023-Jun-15: set albedo to NaN when sza >= 90
#
#######################################################################################################################################################################################################
"""

    canopy_radiation!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Updates canopy radiation profiles for shortwave and longwave radiation, given
- `config` Configurations of spac model
- `spac` `MultiLayerSPAC` type SPAC

"""
canopy_radiation!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; ANGLES, CANOPY, LEAVES, METEO, SOIL) = spac;
    (; DIM_LAYER) = config;

    soil_albedo!(config, SOIL);
    # TODO: note here that this will disable the optical properties of longwave radiation and result in bugs
    if ANGLES.sza < 89
        canopy_optical_properties!(config, CANOPY, ANGLES);
        canopy_optical_properties!(config, CANOPY, LEAVES, SOIL);
        shortwave_radiation!(config, CANOPY, LEAVES, METEO.rad_sw, SOIL);
    else
        CANOPY.RADIATION.r_net_sw .= 0;
        SOIL.ALBEDO.r_net_sw = 0;
        CANOPY.RADIATION.par_in_diffuse = 0;
        CANOPY.RADIATION.par_in_direct = 0;
        CANOPY.RADIATION.par_in = 0;
        CANOPY.RADIATION.par_shaded .= 0;
        CANOPY.RADIATION.par_sunlit .= 0;
        CANOPY.RADIATION.apar_shaded .= 0;
        CANOPY.RADIATION.apar_sunlit .= 0;
        CANOPY.RADIATION.e_v .= 0;
        CANOPY.RADIATION.e_o .= 0;
        CANOPY.RADIATION.albedo .= NaN;

        for _i in 1:DIM_LAYER
            # PPAR for leaves
            LEAVES[_i].ppar_shaded = 0;
            LEAVES[_i].ppar_sunlit .= 0;
        end;
    end;
    longwave_radiation!(CANOPY, LEAVES, METEO.rad_lw, SOIL);

    return nothing
);
