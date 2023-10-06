# Soil albedo values from CLM
const SOIL_ALBEDOS = [0.36 0.61 0.25 0.50;
                      0.34 0.57 0.23 0.46;
                      0.32 0.53 0.21 0.42;
                      0.31 0.51 0.20 0.40;
                      0.30 0.49 0.19 0.38;
                      0.29 0.48 0.18 0.36;
                      0.28 0.45 0.17 0.34;
                      0.27 0.43 0.16 0.32;
                      0.26 0.41 0.15 0.30;
                      0.25 0.39 0.14 0.28;
                      0.24 0.37 0.13 0.26;
                      0.23 0.35 0.12 0.24;
                      0.22 0.33 0.11 0.22;
                      0.20 0.31 0.10 0.20;
                      0.18 0.29 0.09 0.18;
                      0.16 0.27 0.08 0.16;
                      0.14 0.25 0.07 0.14;
                      0.12 0.23 0.06 0.12;
                      0.10 0.21 0.05 0.10;
                      0.08 0.16 0.04 0.08];


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-14: migrate the function from CanopyLayers
#     2022-Jun-14: add method to update broadband or hyperspectral soil albedo
#     2022-Jul-27: add albedo._θ control to HyperspectralSoilAlbedo method (fitting required)
#     2023-Jun-15: make albedo._θ control judge to 0.001
#
#######################################################################################################################################################################################################
"""

    soil_albedo!(config::SPACConfiguration{FT}, sbulk::SoilBulk{FT}, soil::SoilLayer{FT}) where {FT}

Updates lower soil boundary reflectance, given
- `config` Configurations of spac model
- `sbulk` `SoilBulk` type struct
- `soil` `SoilLayer` type struct

"""
function soil_albedo!(config::SPACConfiguration{FT}, sbulk::SoilBulk{FT}, soil::SoilLayer{FT}) where {FT}
    (; SPECTRA, α_CLM, α_FITTING) = config;
    @assert 1 <= sbulk.state.color <=20;

    # if the change of swc is lower than 0.001, do nothing
    if abs(soil.state.θ - sbulk.auxil._θ) < 0.001
        return nothing
    end;

    # use CLM method or Yujie's method
    _rwc::FT = soil.state.θ / soil.state.vc.Θ_SAT;
    _par::FT = SOIL_ALBEDOS[sbulk.state.color,1] * (1 - _rwc) + _rwc * SOIL_ALBEDOS[sbulk.state.color,3];
    _nir::FT = SOIL_ALBEDOS[sbulk.state.color,2] * (1 - _rwc) + _rwc * SOIL_ALBEDOS[sbulk.state.color,4];

    if α_CLM
        _delta = max(0, FT(0.11) - FT(0.4) * soil.state.θ);
        _par = max(SOIL_ALBEDOS[sbulk.state.color,1], SOIL_ALBEDOS[sbulk.state.color,3] + _delta);
        _nir = max(SOIL_ALBEDOS[sbulk.state.color,2], SOIL_ALBEDOS[sbulk.state.color,4] + _delta);
    end;

    # if fitting is disabled, use broadband directly
    if !α_FITTING
        sbulk.auxil.ρ_sw[SPECTRA.IΛ_PAR] .= _par;
        sbulk.auxil.ρ_sw[SPECTRA.IΛ_NIR] .= _nir;

        return nothing
    end;

    # make an initial guess of the weights
    _ρ_sw = similar(sbulk.auxil.ρ_sw);
    _ρ_sw[SPECTRA.IΛ_PAR] .= _par;
    _ρ_sw[SPECTRA.IΛ_NIR] .= _nir;
    sbulk.auxil.weight .= pinv(SPECTRA.MAT_SOIL) * _ρ_sw;

    # function to solve for weights
    @inline _fit(x::Vector{FT}) where {FT} = (
        mul!(_ρ_sw, SPECTRA.MAT_SOIL, x);
        _tmp_vec_nir = abs.(view(_ρ_sw,SPECTRA.IΛ_NIR) .- _nir);
        _diff = ( mean( view(_ρ_sw,SPECTRA.IΛ_PAR) ) - _par ) ^ 2 + mean( _tmp_vec_nir ) ^ 2;

        return -_diff
    );

    # solve for weights
    _ms  = ReduceStepMethodND{FT}(x_mins = FT[-2,-2,-2,-2], x_maxs = FT[2,2,2,2], x_inis = sbulk.auxil.weight, Δ_inis = FT[0.1,0.1,0.1,0.1]);
    _tol = SolutionToleranceND{FT}(FT[0.001,0.001,0.001,0.001], 50);
    _sol = find_peak(_fit, _ms, _tol);
    sbulk.auxil.weight .= _sol;

    # update vectors in soil
    mul!(sbulk.auxil.ρ_sw, SPECTRA.MAT_SOIL, sbulk.auxil.weight);

    # update the albedo._θ to avoid calling this function too many times
    sbulk.auxil._θ = soil.state.θ;

    return nothing
end;
