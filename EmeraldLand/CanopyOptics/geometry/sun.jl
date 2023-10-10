# This file contains function to compute the sun geometry related parameters (apart from sensor geometry)

function sun_geometry!(can::MultiLayerCanopy{FT}) where {FT}
    (; DIM_LAYER, Θ_AZI, _1_AZI, _COS_Θ_AZI, _COS²_Θ_INCL, _COS²_Θ_INCL_AZI) = config;
    (; HOT_SPOT, OPTICS, P_INCL) = can;

    if can.lai == 0
        return nothing
    end;

    # 1. update the canopy optical properties related to extinction and scattering coefficients
    extinction_scattering_coefficients!(config, can);

    can.sensor_geometry.auxil.ko = P_INCL' * OPTICS._ko;
    can.sun_geometry.auxil.ks = P_INCL' * OPTICS._ks;
    OPTICS._bf = P_INCL' * _COS²_Θ_INCL;

    can.sun_geometry.auxil.ddb = (1 + OPTICS._bf) / 2;
    can.sun_geometry.auxil.ddf = (1 - OPTICS._bf) / 2;
    can.sun_geometry.auxil.sdb = (can.sun_geometry.auxil.ks + OPTICS._bf) / 2;
    can.sun_geometry.auxil.sdf = (can.sun_geometry.auxil.ks - OPTICS._bf) / 2;
    can.sensor_geometry.auxil.dob = (can.sensor_geometry.auxil.ko + OPTICS._bf) / 2;
    can.sensor_geometry.auxil.dof = (can.sensor_geometry.auxil.ko - OPTICS._bf) / 2;
    can.sensor_geometry.auxil.sob = P_INCL' * OPTICS._sb;
    can.sensor_geometry.auxil.sof = P_INCL' * OPTICS._sf;

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
    can.sensor_geometry.auxil.po = exp.(can.sensor_geometry.auxil.ko .* can.ci * can.lai .* can._x_bnds);
    can.sun_geometry.auxil.ps = exp.(can.sun_geometry.auxil.ks .* can.ci * can.lai .* can._x_bnds);
    for i in 1:DIM_LAYER
        _ilai = can.ci * can.δlai[i];
        can.sun_geometry.auxil.p_sunlit[i] = 1 / (can.sun_geometry.auxil.ks * _ilai) *
                                             (exp(can.sun_geometry.auxil.ks * can.ci * can.lai * can._x_bnds[i]) - exp(can.sun_geometry.auxil.ks * can.ci * can.lai * can._x_bnds[i+1]));
    end;

    _dso = sqrt( tand(can.sun_geometry.state.sza) ^ 2 +
                 tand(can.sensor_geometry.state.vza) ^ 2 -
                 2 * tand(can.sun_geometry.state.sza) * tand(can.sensor_geometry.state.vza) * cosd(can.sensor_geometry.state.vaa - can.sun_geometry.state.saa) );
    @inline _pdf(x::FT) where {FT} = (
        _Σk = can.sensor_geometry.auxil.ko + can.sun_geometry.auxil.ks;
        _Πk = can.sensor_geometry.auxil.ko * can.sun_geometry.auxil.ks;
        _cl = can.ci * can.lai;
        _α  = _dso / HOT_SPOT * 2 / _Σk;

        if _dso == 0
            return exp( (_Σk - sqrt(_Πk)) * _cl * x )
        end;

        return exp( _Σk * _cl * x + sqrt(_Πk) * _cl / _α * (1 - exp(_α * x)) )
    );

    for i in eachindex(can._x_bnds)
        can.sensor_geometry.auxil.pso[i] = quadgk(_pdf, can._x_bnds[i] - FT(1)/DIM_LAYER, can._x_bnds[i]; rtol = 1e-2)[1] * DIM_LAYER;
    end;

    return nothing
end;
