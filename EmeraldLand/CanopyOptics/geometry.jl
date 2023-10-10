#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-07: migrate the function from CanopyLayers
#     2022-Jun-09: rename function to canopy_optical_properties!
#
#######################################################################################################################################################################################################
"""
This function updates canopy optical properties for canopy. The supported methods are to
- Update the extinction coefficients
- Update the soil boundary conditions (not public function)
- Update scattering coefficient matrices

"""
function canopy_optical_properties! end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-07: migrate the function from CanopyLayers
#     2022-Jun-09: rename function to canopy_optical_properties!
#     2022-Jun-09: update p_sunlit
#     2023-Mar-11: add code to account for the case of LAI == 0
#     2023-Apr-25: use exp to compute po, ps, and p_sunlit
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

    canopy_optical_properties!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT}

Updates canopy optical properties (extinction coefficients for direct and diffuse light) based on the SAIL model, given
- `config` SPAC configurations
- `can` `MultiLayerCanopy` type struct

"""
canopy_optical_properties!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT} = (
    (; DIM_LAYER, Θ_AZI, _1_AZI, _COS_Θ_AZI, _COS²_Θ_INCL, _COS²_Θ_INCL_AZI) = config;
    (; OPTICS) = can;

    if can.structure.state.lai == 0
        return nothing
    end;

    # 1. update the canopy optical properties related to extinction and scattering coefficients
    extinction_scattering_coefficients!(config, can);

    can.sensor_geometry.auxil.ko = can.structure.state.p_incl' * OPTICS._ko;
    can.sun_geometry.auxil.ks = can.structure.state.p_incl' * OPTICS._ks;
    OPTICS._bf = can.structure.state.p_incl' * _COS²_Θ_INCL;

    can.sun_geometry.auxil.ddb = (1 + OPTICS._bf) / 2;
    can.sun_geometry.auxil.ddf = (1 - OPTICS._bf) / 2;
    can.sun_geometry.auxil.sdb = (can.sun_geometry.auxil.ks + OPTICS._bf) / 2;
    can.sun_geometry.auxil.sdf = (can.sun_geometry.auxil.ks - OPTICS._bf) / 2;
    can.sensor_geometry.auxil.dob = (can.sensor_geometry.auxil.ko + OPTICS._bf) / 2;
    can.sensor_geometry.auxil.dof = (can.sensor_geometry.auxil.ko - OPTICS._bf) / 2;
    can.sensor_geometry.auxil.sob = can.structure.state.p_incl' * OPTICS._sb;
    can.sensor_geometry.auxil.sof = can.structure.state.p_incl' * OPTICS._sf;

    # 2. update the matrices fs and fo
    OPTICS._cos_θ_azi_raa .= cosd.(Θ_AZI .- (can.sensor_geometry.state.vaa - can.sun_geometry.state.saa));
    mul!(OPTICS._tmp_mat_incl_azi_1, OPTICS._Co, _1_AZI');
    mul!(OPTICS._tmp_mat_incl_azi_2, OPTICS._So, OPTICS._cos_θ_azi_raa');
    OPTICS.fo .= (OPTICS._tmp_mat_incl_azi_1 .+ OPTICS._tmp_mat_incl_azi_2) ./ cosd(can.sensor_geometry.state.vza);
    OPTICS._abs_fo .= abs.(OPTICS.fo);

    mul!(OPTICS._tmp_mat_incl_azi_1, OPTICS._Cs, _1_AZI');
    mul!(OPTICS._tmp_mat_incl_azi_2, OPTICS._Ss, _COS_Θ_AZI');
    OPTICS.fs .= (OPTICS._tmp_mat_incl_azi_1 .+ OPTICS._tmp_mat_incl_azi_2) ./ cosd(can.sun_geometry.state.sza);
    OPTICS._abs_fs .= abs.(OPTICS.fs);

    OPTICS._fo_cos_θ_incl .= OPTICS.fo .* _COS²_Θ_INCL_AZI;
    OPTICS._fs_cos_θ_incl .= OPTICS.fs .* _COS²_Θ_INCL_AZI;
    OPTICS._fs_fo .= OPTICS.fs .* OPTICS.fo;
    OPTICS._abs_fs_fo .= abs.(OPTICS._fs_fo);

    # 3. update the viewing fraction ps, po, pso, and p_sunlit
    can.sensor_geometry.auxil.po = exp.(can.sensor_geometry.auxil.ko .* can.structure.auxil.ci * can.structure.state.lai .* can.structure.auxil.x_bnds);
    can.sun_geometry.auxil.ps = exp.(can.sun_geometry.auxil.ks .* can.structure.auxil.ci * can.structure.state.lai .* can.structure.auxil.x_bnds);
    for i in 1:DIM_LAYER
        _ilai = can.structure.auxil.ci * can.structure.state.δlai[i];
        can.sun_geometry.auxil.p_sunlit[i] = 1 / (can.sun_geometry.auxil.ks * _ilai) *
                                             (exp(can.sun_geometry.auxil.ks * can.structure.auxil.ci * can.structure.state.lai * can.structure.auxil.x_bnds[i]) -
                                              exp(can.sun_geometry.auxil.ks * can.structure.auxil.ci * can.structure.state.lai * can.structure.auxil.x_bnds[i+1]));
    end;

    _dso = sqrt( tand(can.sun_geometry.state.sza) ^ 2 +
                 tand(can.sensor_geometry.state.vza) ^ 2 -
                 2 * tand(can.sun_geometry.state.sza) * tand(can.sensor_geometry.state.vza) * cosd(can.sensor_geometry.state.vaa - can.sun_geometry.state.saa) );
    @inline _pdf(x::FT) where {FT} = (
        _Σk = can.sensor_geometry.auxil.ko + can.sun_geometry.auxil.ks;
        _Πk = can.sensor_geometry.auxil.ko * can.sun_geometry.auxil.ks;
        _cl = can.structure.auxil.ci * can.structure.state.lai;
        _α  = _dso / can.structure.state.hot_spot * 2 / _Σk;

        if _dso == 0
            return exp( (_Σk - sqrt(_Πk)) * _cl * x )
        end;

        return exp( _Σk * _cl * x + sqrt(_Πk) * _cl / _α * (1 - exp(_α * x)) )
    );

    for i in eachindex(can.structure.auxil.x_bnds)
        can.sensor_geometry.auxil.pso[i] = quadgk(_pdf, can.structure.auxil.x_bnds[i] - FT(1)/DIM_LAYER, can.structure.auxil.x_bnds[i]; rtol = 1e-2)[1] * DIM_LAYER;
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-08: migrate the function from CanopyLayers
#     2022-Jun-09: rename the function from canopy_matrices! to canopy_optical_properties!
#     2022-Jun-09: move part of the short_wave! code into canopy_optical_properties!
#     2022-Jun-10: add function to compute longwave reflectance, transmittance, and emissivity
#     2022-Jun-29: use Leaf for the hyperspectral RT
#     2023-Mar-11: add code to account for the case of LAI == 0
#     2023-Apr-25: use exp to computethe scattering coefficients
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

    canopy_optical_properties!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}, leaves::Vector{Leaf{FT}}, sbulk::SoilBulk{FT}) where {FT}

Updates canopy optical properties (scattering coefficient matrices), given
- `config` Configuration for `MultiLayerSPAC`
- `can` `MultiLayerCanopy` type struct
- `leaves` Vector of `Leaf`
- `sbulk` `SoilBulk` type struct

"""
canopy_optical_properties!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}, leaves::Vector{Leaf{FT}}, sbulk::SoilBulk{FT}) where {FT} = (
    (; DIM_LAYER) = config;
    (; OPTICS) = can;
    @assert length(leaves) == DIM_LAYER "Number of leaves must be equal to the canopy layers!";

    if can.structure.state.lai == 0
        OPTICS.ρ_dd  .= 0;
        OPTICS.ρ_lw  .= 0;
        OPTICS.ρ_sd  .= 0;
        OPTICS.τ_dd  .= 0;
        OPTICS.τ_lw  .= 0;
        OPTICS.τ_sd  .= 0;
        OPTICS._τ_ss .= 0;
        OPTICS.ρ_dd[:,end] .= sbulk.auxil.ρ_sw;
        OPTICS.ρ_sd[:,end] .= sbulk.auxil.ρ_sw;
        OPTICS.ρ_lw[end] = sbulk.auxil.ρ_lw;

        return nothing
    end;

    # 1. update the scattering coefficients for different layers
    for i in eachindex(leaves)
        OPTICS.σ_ddb[:,i] .= can.sun_geometry.auxil.ddb * leaves[i].bio.auxil.ρ_leaf .+ can.sun_geometry.auxil.ddf * leaves[i].bio.auxil.τ_leaf;
        OPTICS.σ_ddf[:,i] .= can.sun_geometry.auxil.ddf * leaves[i].bio.auxil.ρ_leaf .+ can.sun_geometry.auxil.ddb * leaves[i].bio.auxil.τ_leaf;
        OPTICS.σ_sdb[:,i] .= can.sun_geometry.auxil.sdb * leaves[i].bio.auxil.ρ_leaf .+ can.sun_geometry.auxil.sdf * leaves[i].bio.auxil.τ_leaf;
        OPTICS.σ_sdf[:,i] .= can.sun_geometry.auxil.sdf * leaves[i].bio.auxil.ρ_leaf .+ can.sun_geometry.auxil.sdb * leaves[i].bio.auxil.τ_leaf;
        OPTICS.σ_dob[:,i] .= can.sensor_geometry.auxil.dob * leaves[i].bio.auxil.ρ_leaf .+ can.sensor_geometry.auxil.dof * leaves[i].bio.auxil.τ_leaf;
        OPTICS.σ_dof[:,i] .= can.sensor_geometry.auxil.dof * leaves[i].bio.auxil.ρ_leaf .+ can.sensor_geometry.auxil.dob * leaves[i].bio.auxil.τ_leaf;
        OPTICS.σ_so[:,i]  .= can.sensor_geometry.auxil.sob * leaves[i].bio.auxil.ρ_leaf .+ can.sensor_geometry.auxil.sof * leaves[i].bio.auxil.τ_leaf;
    end;

    # 2. update the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    OPTICS._τ_ss .= exp.(-1 .* can.sun_geometry.auxil.ks .* can.structure.state.δlai .* can.structure.auxil.ci);
    OPTICS._τ_dd .= exp.(-1 .* (1 .- OPTICS.σ_ddf) .* can.structure.state.δlai' .* can.structure.auxil.ci);
    OPTICS._τ_sd .= 1 .- exp.(-1 .* OPTICS.σ_sdf .* can.structure.state.δlai' .* can.structure.auxil.ci);
    OPTICS._ρ_dd .= 1 .- exp.(-1 .* OPTICS.σ_ddb .* can.structure.state.δlai' .* can.structure.auxil.ci);
    OPTICS._ρ_sd .= 1 .- exp.(-1 .* OPTICS.σ_sdb .* can.structure.state.δlai' .* can.structure.auxil.ci);

    # 3. update the effective reflectance per layer
    OPTICS.ρ_dd[:,end] .= sbulk.auxil.ρ_sw;
    OPTICS.ρ_sd[:,end] .= sbulk.auxil.ρ_sw;

    for i in DIM_LAYER:-1:1
        _r_dd__ = view(OPTICS._ρ_dd,:,i  );    # reflectance without correction
        _r_dd_i = view(OPTICS.ρ_dd ,:,i  );    # reflectance of the upper boundary (i)
        _r_dd_j = view(OPTICS.ρ_dd ,:,i+1);    # reflectance of the lower boundary (i+1)
        _r_sd__ = view(OPTICS._ρ_sd,:,i  );    # reflectance without correction
        _r_sd_i = view(OPTICS.ρ_sd ,:,i  );    # reflectance of the upper boundary (i)
        _r_sd_j = view(OPTICS.ρ_sd ,:,i+1);    # reflectance of the lower boundary (i+1)
        _t_dd__ = view(OPTICS._τ_dd,:,i  );    # transmittance without correction
        _t_dd_i = view(OPTICS.τ_dd ,:,i  );    # transmittance of the layer (i)
        _t_sd__ = view(OPTICS._τ_sd,:,i  );    # transmittance without correction
        _t_sd_i = view(OPTICS.τ_sd ,:,i  );    # transmittance of the layer (i)
        _t_ss__ = view(OPTICS._τ_ss,i);        # transmittance for directional->directional

        OPTICS._tmp_vec_λ .= 1 .- _r_dd__ .* _r_dd_j;

        _t_sd_i .= (_t_sd__ .+ _t_ss__ .* _r_sd_j .* _r_dd__) ./ OPTICS._tmp_vec_λ;             # it + it-jr-ir, then rescale
        _t_dd_i .= _t_dd__ ./ OPTICS._tmp_vec_λ;                                                # it, then rescale
        _r_sd_i .= _r_sd__ .+ _t_ss__ .* _r_sd_j .* _t_dd__ .+ _t_sd_i .* _r_dd_j .* _t_dd__;   # ir + it-jr-it(v) + it-jr_dd-it
        _r_dd_i .= _r_dd__ .+ _t_dd__ .* _r_dd_j .* _t_dd_i;                                    # ir + it-jr-it
    end;

    # 4. compute longwave effective emissivity, reflectance, and transmittance (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    for i in 1:DIM_LAYER
        _ilai = can.structure.state.δlai[i] * can.structure.auxil.ci;
        _σ_lw_b = can.sun_geometry.auxil.ddb * leaves[i].bio.auxil.ρ_LW + can.sun_geometry.auxil.ddf * leaves[i].bio.auxil.τ_LW;
        _σ_lw_f = can.sun_geometry.auxil.ddf * leaves[i].bio.auxil.ρ_LW + can.sun_geometry.auxil.ddb * leaves[i].bio.auxil.τ_LW;
        OPTICS._τ_lw[i] = exp(-1 * (1 - _σ_lw_f) * _ilai);
        OPTICS._ρ_lw[i] = 1 - exp(-1 * _σ_lw_b * _ilai);
        OPTICS.ϵ[i] = 1 - OPTICS._τ_lw[i] - OPTICS._ρ_lw[i];
    end;

    # 5. update the effective longwave reflectance and transmittance
    OPTICS.ρ_lw[end] = sbulk.auxil.ρ_lw;

    for i in DIM_LAYER:-1:1
        _dnorm = 1 - OPTICS._ρ_lw[i] * OPTICS.ρ_lw[i+1];

        OPTICS.τ_lw[i] = OPTICS._τ_lw[i] / _dnorm;                                                    # it, then rescale
        OPTICS.ρ_lw[i] = OPTICS._ρ_lw[i] + OPTICS._τ_lw[i] * OPTICS.ρ_lw[i+1] * OPTICS.τ_lw[i];    # ir + it-jr-it
    end;

    return nothing
);
