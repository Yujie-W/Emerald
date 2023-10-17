# This file contains function to compute soil albedo

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

    soil_albedo!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Updates lower soil boundary reflectance, given
- `config` Configurations of spac model
- `spac` SPAC

"""
function soil_albedo!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    (; SPECTRA, α_CLM, α_FITTING) = config;
    (; SOIL_BULK, SOILS) = spac;

    @assert 1 <= SOIL_BULK.state.color <=20;

    # if the change of swc is lower than 0.001, do nothing
    if abs(SOILS[1].state.θ - SOIL_BULK.auxil._θ) < 0.001
        return nothing
    end;

    # use linear interpolation method or CLM method (with upper limit)
    rwc = SOILS[1].state.θ / SOILS[1].state.vc.Θ_SAT;
    par::FT = SOIL_ALBEDOS[SOIL_BULK.state.color,1] * (1 - rwc) + rwc * SOIL_ALBEDOS[SOIL_BULK.state.color,3];
    nir::FT = SOIL_ALBEDOS[SOIL_BULK.state.color,2] * (1 - rwc) + rwc * SOIL_ALBEDOS[SOIL_BULK.state.color,4];

    if α_CLM
        delta = max(0, FT(0.11) - FT(0.4) * SOILS[1].state.θ);
        par = max(SOIL_ALBEDOS[SOIL_BULK.state.color,1], SOIL_ALBEDOS[SOIL_BULK.state.color,3] + delta);
        nir = max(SOIL_ALBEDOS[SOIL_BULK.state.color,2], SOIL_ALBEDOS[SOIL_BULK.state.color,4] + delta);
    end;

    # if fitting is disabled, use broadband directly
    if !α_FITTING
        SOIL_BULK.auxil.ρ_sw[SPECTRA.IΛ_PAR] .= par;
        SOIL_BULK.auxil.ρ_sw[SPECTRA.IΛ_NIR] .= nir;

        return nothing
    end;

    #
    # TODO: use a new soil moddel for this, do not GSV which is not process-based
    #
    # make an initial guess of the weights
    ρ_sw = similar(SOIL_BULK.auxil.ρ_sw);
    ρ_sw[SPECTRA.IΛ_PAR] .= par;
    ρ_sw[SPECTRA.IΛ_NIR] .= nir;
    SOIL_BULK.auxil.weight .= pinv(SPECTRA.MAT_SOIL) * ρ_sw;

    # function to solve for weights
    @inline _fit(x::Vector{FT}) where {FT} = (
        mul!(ρ_sw, SPECTRA.MAT_SOIL, x);
        tmp_vec_nir = abs.(view(ρ_sw,SPECTRA.IΛ_NIR) .- nir);
        diff = ( mean( view(ρ_sw,SPECTRA.IΛ_PAR) ) - par ) ^ 2 + mean( tmp_vec_nir ) ^ 2;

        return -diff
    );

    # solve for weights
    ms = ReduceStepMethodND{FT}(x_mins = FT[-2,-2,-2,-2], x_maxs = FT[2,2,2,2], x_inis = SOIL_BULK.auxil.weight, Δ_inis = FT[0.1,0.1,0.1,0.1]);
    tol = SolutionToleranceND{FT}(FT[0.001,0.001,0.001,0.001], 50);
    sol = find_peak(_fit, ms, tol);
    SOIL_BULK.auxil.weight .= sol;

    # update vectors in soil
    mul!(SOIL_BULK.auxil.ρ_sw, SPECTRA.MAT_SOIL, SOIL_BULK.auxil.weight);

    # update the albedo._θ to avoid calling this function too many times
    SOIL_BULK.auxil._θ = SOILS[1].state.θ;

    return nothing
end;
