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
function shortwave_radiation! end;


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

    shortwave_radiation!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::ShortwaveRadiation{FT}, sbulk::SoilBulk{FT}) where {FT}

Updates canopy radiation profiles for shortwave radiation, given
- `config` Configuration for `MultiLayerSPAC`
- `can` `MultiLayerCanopy` type struct
- `leaves` Vector of `Leaf`
- `rad` Incoming shortwave radiation
- `sbulk` Soil bulk parameters

"""
shortwave_radiation!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::ShortwaveRadiation{FT}, sbulk::SoilBulk{FT}) where {FT} = (
    (; DIM_LAYER, SPECTRA) = config;
    (; OPTICS, P_INCL, RADIATION) = can;

    if can.lai == 0
        for i in 1:(DIM_LAYER+1)
            RADIATION.e_direct[:,i] .= rad.e_direct;
            RADIATION.e_diffuse_down[:,i] .= rad.e_diffuse;
            RADIATION.e_diffuse_up[:,i] .= 0;
            RADIATION.e_v[:,i] .= 0;
        end;
        RADIATION.e_diffuse_up[:,end;] .= view(OPTICS.ρ_sd,:,DIM_LAYER+1) .* view(RADIATION.e_direct,:,DIM_LAYER+1) .+ view(OPTICS.ρ_dd,:,DIM_LAYER+1) .* view(RADIATION.e_diffuse_down,:,DIM_LAYER+1);
        RADIATION.e_v[:,end;] .= view(RADIATION.e_diffuse_up,:,DIM_LAYER+1);
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

    for i in 1:DIM_LAYER
        _e_d_i = view(RADIATION.e_diffuse_down,:,i  );     # downward diffuse radiation at upper boundary
        _e_d_j = view(RADIATION.e_diffuse_down,:,i+1);     # downward diffuse radiation at lower boundary
        _e_s_i = view(RADIATION.e_direct      ,:,i  );     # direct radiation at upper boundary
        _e_s_j = view(RADIATION.e_direct      ,:,i+1);     # direct radiation at lower boundary
        _e_u_i = view(RADIATION.e_diffuse_up  ,:,i  );     # upward diffuse radiation at upper boundary

        _r_dd_i = view(OPTICS.ρ_dd ,:,i);  # reflectance of the upper boundary (i)
        _r_sd_i = view(OPTICS.ρ_sd ,:,i);  # reflectance of the upper boundary (i)
        _t_dd_i = view(OPTICS.τ_dd ,:,i);  # transmittance of the layer (i)
        _t_sd_i = view(OPTICS.τ_sd ,:,i);  # transmittance of the layer (i)
        _t_ss__ = view(OPTICS._τ_ss,i);    # transmittance for directional->directional

        _e_s_j .= _t_ss__ .* _e_s_i;
        _e_d_j .= _t_sd_i .* _e_s_i .+ _t_dd_i .* _e_d_i;
        _e_u_i .= _r_sd_i .* _e_s_i .+ _r_dd_i .* _e_d_i;
    end;

    RADIATION.e_diffuse_up[:,end;] .= view(OPTICS.ρ_sd,:,DIM_LAYER+1) .* view(RADIATION.e_direct,:,DIM_LAYER+1) .+ view(OPTICS.ρ_dd,:,DIM_LAYER+1) .* view(RADIATION.e_diffuse_down,:,DIM_LAYER+1);

    # 2. update the sunlit and shaded sum radiation and total absorbed radiation per layer and for soil
    for i in 1:DIM_LAYER
        _a_s_i = view(RADIATION.e_net_direct  ,:,i  );     # net absorbed direct radiation
        _a_d_i = view(RADIATION.e_net_diffuse ,:,i  );     # net absorbed diffuse radiation
        _e_d_i = view(RADIATION.e_diffuse_down,:,i  );     # downward diffuse radiation at upper boundary
        _e_s_i = view(RADIATION.e_direct      ,:,i  );     # direct radiation at upper boundary
        _e_u_j = view(RADIATION.e_diffuse_up  ,:,i+1);     # upward diffuse radiation at lower boundary
        _p_s_i = view(RADIATION.e_sum_direct  ,:,i  );     # sum direct radiation
        _p_d_i = view(RADIATION.e_sum_diffuse ,:,i  );     # sum diffuse radiation

        _r_dd__ = view(OPTICS._ρ_dd,:,i);  # reflectance of the upper boundary (i)
        _r_sd__ = view(OPTICS._ρ_sd,:,i);  # reflectance of the upper boundary (i)
        _t_dd__ = view(OPTICS._τ_dd,:,i);  # transmittance of the layer (i)
        _t_sd__ = view(OPTICS._τ_sd,:,i);  # transmittance of the layer (i)
        _t_ss__ = view(OPTICS._τ_ss,i);    # transmittance for directional->directional

        _p_s_i .= _e_s_i;
        _p_d_i .= _e_d_i .+ _e_u_j;

        _a_s_i .= _p_s_i .* (1 .- _t_ss__ .- _t_sd__ .- _r_sd__);
        _a_d_i .= _p_d_i .* (1 .- _t_dd__ .- _r_dd__);
    end;

    # 3. compute the spectra at the observer direction
    for i in 1:DIM_LAYER
        _e_d_i = view(RADIATION.e_diffuse_down,:,i);   # downward diffuse radiation at upper boundary
        _e_u_i = view(RADIATION.e_diffuse_up  ,:,i);   # upward diffuse radiation at upper boundary
        _e_v_i = view(RADIATION.e_v,:,i);

        _dob_i = view(OPTICS.σ_dob,:,i);   # scattering coefficient backward for diffuse->observer
        _dof_i = view(OPTICS.σ_dof,:,i);   # scattering coefficient forward for diffuse->observer
        _so__i = view(OPTICS.σ_so ,:,i);   # bidirectional from solar to observer

        _e_v_i  .= can.sensor_geometry.auxil.po[i] .* _dob_i .* _e_d_i;
        _e_v_i .+= can.sensor_geometry.auxil.po[i] .* _dof_i .* _e_u_i;
        _e_v_i .+= can.sensor_geometry.auxil.pso[i] .* _so__i .* rad.e_direct;
        _e_v_i .*= can.δlai[i] * can.ci;
    end;
    RADIATION.e_v[:,end;] .= can.sensor_geometry.auxil.po[end;] .* view(RADIATION.e_diffuse_up,:,DIM_LAYER+1);

    for i in eachindex(RADIATION.e_o)
        RADIATION.e_o[i] = sum(view(RADIATION.e_v,i,:)) / FT(pi);
    end;

    RADIATION.albedo .= RADIATION.e_o * FT(pi) ./ (rad.e_direct .+ rad.e_diffuse);

    # 4. compute net absorption for leaves and soil
    for i in 1:DIM_LAYER
        _Σ_shaded = view(RADIATION.e_net_diffuse,:,i)' * SPECTRA.ΔΛ / 1000 / can.δlai[i];
        _Σ_sunlit = view(RADIATION.e_net_direct ,:,i)' * SPECTRA.ΔΛ / 1000 / can.δlai[i];
        RADIATION.r_net_sw_shaded[i] = _Σ_shaded;
        RADIATION.r_net_sw_sunlit[i] = _Σ_sunlit / can.sun_geometry.auxil.p_sunlit[i] + _Σ_shaded;
        RADIATION.r_net_sw[i] = _Σ_shaded * (1 - can.sun_geometry.auxil.p_sunlit[i]) + _Σ_sunlit * can.sun_geometry.auxil.p_sunlit[i];
    end;

    sbulk.auxil.e_net_direct .= view(RADIATION.e_direct,:,DIM_LAYER+1) .* (1 .- sbulk.auxil.ρ_sw);
    sbulk.auxil.e_net_diffuse .= view(RADIATION.e_diffuse_down,:,DIM_LAYER+1) .* (1 .- sbulk.auxil.ρ_sw);
    sbulk.auxil.r_net_sw = (sbulk.auxil.e_net_direct' * SPECTRA.ΔΛ + sbulk.auxil.e_net_diffuse' * SPECTRA.ΔΛ) / 1000;

    # 5. compute top-of-canopy and leaf level PAR, APAR, and PPAR per ground area
    RADIATION._par_shaded .= photon.(SPECTRA.Λ_PAR, view(rad.e_diffuse,SPECTRA.IΛ_PAR)) .* 1000;
    RADIATION._par_sunlit .= photon.(SPECTRA.Λ_PAR, view(rad.e_direct ,SPECTRA.IΛ_PAR)) .* 1000;
    RADIATION.par_in_diffuse = RADIATION._par_shaded' * SPECTRA.ΔΛ_PAR;
    RADIATION.par_in_direct = RADIATION._par_sunlit' * SPECTRA.ΔΛ_PAR;
    RADIATION.par_in = RADIATION.par_in_diffuse + RADIATION.par_in_direct;

    mul!(OPTICS._tmp_vec_azi, OPTICS._abs_fs', P_INCL);
    _normi = 1 / mean(OPTICS._tmp_vec_azi);

    for i in 1:DIM_LAYER
        _α_apar = view(leaves[i].bio.auxil.f_ppar, SPECTRA.IΛ_PAR);

        # convert energy to quantum unit for PAR, APAR and PPAR per leaf area
        RADIATION._par_shaded  .= photon.(SPECTRA.Λ_PAR, view(RADIATION.e_sum_diffuse,SPECTRA.IΛ_PAR,i)) .* 1000;
        RADIATION._par_sunlit  .= photon.(SPECTRA.Λ_PAR, view(RADIATION.e_sum_direct ,SPECTRA.IΛ_PAR,i)) .* 1000 ./ can.sun_geometry.auxil.p_sunlit[i];
        RADIATION._apar_shaded .= photon.(SPECTRA.Λ_PAR, view(RADIATION.e_net_diffuse,SPECTRA.IΛ_PAR,i)) .* 1000 ./ can.δlai[i];
        RADIATION._apar_sunlit .= photon.(SPECTRA.Λ_PAR, view(RADIATION.e_net_direct ,SPECTRA.IΛ_PAR,i)) .* 1000 ./ can.δlai[i] ./ can.sun_geometry.auxil.p_sunlit[i];
        RADIATION._ppar_shaded .= RADIATION._apar_shaded .* _α_apar;
        RADIATION._ppar_sunlit .= RADIATION._apar_sunlit .* _α_apar;

        # PAR for leaves
        _Σ_par_dif = RADIATION._par_shaded' * SPECTRA.ΔΛ_PAR;
        _Σ_par_dir = RADIATION._par_sunlit' * SPECTRA.ΔΛ_PAR * _normi;
        RADIATION.par_shaded[i] = _Σ_par_dif;
        RADIATION.par_sunlit[:,:,i] .= OPTICS._abs_fs_fo .* _Σ_par_dir;
        RADIATION.par_sunlit[:,:,i] .+= _Σ_par_dif;

        # APAR for leaves
        _Σ_apar_dif = RADIATION._apar_shaded' * SPECTRA.ΔΛ_PAR;
        _Σ_apar_dir = RADIATION._apar_sunlit' * SPECTRA.ΔΛ_PAR * _normi;
        RADIATION.apar_shaded[i] = _Σ_apar_dif;
        RADIATION.apar_sunlit[:,:,i] .= OPTICS._abs_fs_fo .* _Σ_apar_dir;
        RADIATION.apar_sunlit[:,:,i] .+= _Σ_apar_dif;

        # PPAR for leaves
        _Σ_ppar_dif = RADIATION._ppar_shaded' * SPECTRA.ΔΛ_PAR;
        _Σ_ppar_dir = RADIATION._ppar_sunlit' * SPECTRA.ΔΛ_PAR * _normi;
        leaves[DIM_LAYER+1-i].flux.auxil.ppar_shaded  = _Σ_ppar_dif;
        leaves[DIM_LAYER+1-i].flux.auxil.ppar_sunlit .= OPTICS._abs_fs_fo .* _Σ_ppar_dir .+ _Σ_ppar_dif;
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
function longwave_radiation! end;


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

    longwave_radiation!(can::MultiLayerCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::FT, sbulk::SoilBulk{FT}, soil::SoilLayer{FT}) where {FT}

Updates canopy radiation profiles for shortwave or longwave radiation, given
- `can` `MultiLayerCanopy` type struct
- `leaves` Vector of `Leaf`
- `rad` Incoming longwave radiation
- `sbulk` Soil bulk parameters
- `soil` First soil layer

"""
longwave_radiation!(can::MultiLayerCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::FT, sbulk::SoilBulk{FT}, soil::SoilLayer{FT}) where {FT} = (
    (; OPTICS, RADIATION) = can;

    _nlayers = length(leaves);

    if can.lai == 0
        _r_lw_soil = K_STEFAN(FT) * (1 - sbulk.auxil.ρ_lw) * soil.auxil.t ^ 4;
        RADIATION.r_lw .= 0;
        RADIATION.r_net_lw .= 0;
        RADIATION.r_lw_up .= rad * sbulk.auxil.ρ_lw + _r_lw_soil;
        sbulk.auxil.r_net_lw = rad * (1 - sbulk.auxil.ρ_lw) - _r_lw_soil;

        return nothing
    end;

    # 1. compute longwave radiation out from the leaves and soil
    for i in eachindex(leaves)
        RADIATION.r_lw[i] = K_STEFAN(FT) * OPTICS.ϵ[i] * leaves[i].energy.auxil.t ^ 4;
    end;

    _r_lw_soil = K_STEFAN(FT) * (1 - sbulk.auxil.ρ_lw) * soil.auxil.t ^ 4;

    # 2. account for the longwave emission from bottom to up
    RADIATION._r_emit_up[end;] = _r_lw_soil;

    for i in _nlayers:-1:1
        _r__ = OPTICS._ρ_lw[i];
        _r_j = OPTICS.ρ_lw[i+1];
        _t__ = OPTICS._τ_lw[i];

        _dnorm = 1 - _r__ * _r_j;

        RADIATION._r_emit_down[i] = (RADIATION._r_emit_up[i+1] * _r__ + RADIATION.r_lw[i]) / _dnorm;
        RADIATION._r_emit_up[i] = RADIATION._r_emit_down[i] * _r_j * _t__ + RADIATION._r_emit_up[i+1] * _t__ + RADIATION.r_lw[i];
    end;

    # 3. account for the longwave emission from up to bottom
    RADIATION.r_lw_down[1] = rad;

    for i in 1:_nlayers
        _r_i = OPTICS.ρ_lw[i];
        _t_i = OPTICS.τ_lw[i];

        RADIATION.r_lw_down[i+1] = RADIATION.r_lw_down[i] * _t_i + RADIATION._r_emit_down[i];
        RADIATION.r_lw_up[i+1] = RADIATION.r_lw_down[i] * _r_i + RADIATION._r_emit_up[i];
    end;

    RADIATION.r_lw_up[end;] = RADIATION.r_lw_down[end;] * sbulk.auxil.ρ_lw + _r_lw_soil;

    # 4. compute the net longwave radiation per canopy layer and soil
    for i in 1:_nlayers
        RADIATION.r_net_lw[i] = (RADIATION.r_lw_down[i] + RADIATION.r_lw_up[i+1]) * (1 - OPTICS._ρ_lw[i] - OPTICS._τ_lw[i]) - 2 * RADIATION.r_lw[i];
    end;

    sbulk.auxil.r_net_lw = RADIATION.r_lw_down[end;] * (1 - sbulk.auxil.ρ_lw) - _r_lw_soil;

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
function canopy_radiation! end;


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
    (; CANOPY, LEAVES, METEO, SOIL_BULK, SOILS) = spac;
    (; DIM_LAYER) = config;

    soil_albedo!(config, SOIL_BULK, SOILS[1]);
    # TODO: note here that this will disable the optical properties of longwave radiation and result in bugs
    if CANOPY.sun_geometry.state.sza < 89
        canopy_optical_properties!(config, CANOPY);
        canopy_optical_properties!(config, CANOPY, LEAVES, SOIL_BULK);
        shortwave_radiation!(config, CANOPY, LEAVES, METEO.rad_sw, SOIL_BULK);
    else
        CANOPY.RADIATION.r_net_sw .= 0;
        SOIL_BULK.auxil.r_net_sw = 0;
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

        for i in 1:DIM_LAYER
            # PPAR for leaves
            LEAVES[i].flux.auxil.ppar_shaded = 0;
            LEAVES[i].flux.auxil.ppar_sunlit .= 0;
        end;
    end;
    longwave_radiation!(CANOPY, LEAVES, METEO.rad_lw, SOIL_BULK, SOILS[1]);

    return nothing
);
