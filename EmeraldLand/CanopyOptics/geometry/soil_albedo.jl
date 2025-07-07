# This file contains function to compute soil albedo

# Soil albedo values from CLM
#                     PAR  NIR  PAR  NIR
#                     DRY  DRY  WET  WET
const SOIL_ALBEDOS = [0.36 0.61 0.25 0.50;    # color = 1
                      0.34 0.57 0.23 0.46;    # color = 2
                      0.32 0.53 0.21 0.42;    # color = 3
                      0.31 0.51 0.20 0.40;    # color = 4
                      0.30 0.49 0.19 0.38;    # color = 5
                      0.29 0.48 0.18 0.36;    # color = 6
                      0.28 0.45 0.17 0.34;    # color = 7
                      0.27 0.43 0.16 0.32;    # color = 8
                      0.26 0.41 0.15 0.30;    # color = 9
                      0.25 0.39 0.14 0.28;    # color = 10
                      0.24 0.37 0.13 0.26;    # color = 11
                      0.23 0.35 0.12 0.24;    # color = 12
                      0.22 0.33 0.11 0.22;    # color = 13
                      0.20 0.31 0.10 0.20;    # color = 14
                      0.18 0.29 0.09 0.18;    # color = 15
                      0.16 0.27 0.08 0.16;    # color = 16
                      0.14 0.25 0.07 0.14;    # color = 17
                      0.12 0.23 0.06 0.12;    # color = 18
                      0.10 0.21 0.05 0.10;    # color = 19
                      0.08 0.16 0.04 0.08];   # color = 20


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-14: migrate the function from CanopyLayers
#     2022-Jun-14: add method to update broadband or hyperspectral soil albedo
#     2023-Oct-26: add methods for four soil albedo algorithms
#     2024-Apr-19: add new method to prescribe soil albedo (do nothing)
#     2024-Nov-18: add new method to prescribe soil albedo with broadband values and hyperspectral flag
#     2025-Jun-05: account for ice volume in the soil albedo calculation
#     2025-Jun-06: make sure soil albedo is at least 0.01 when using hyperspectral method
# Bug fixes:
#     2025-Jun-06: fix the soil albedo calculation for CLM method (used the typo max instead of min)
#
#######################################################################################################################################################################################################
"""

    soil_albedo!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Updates lower soil boundary reflectance, given
- `config` Configurations of spac model
- `spac` SPAC

"""
function soil_albedo! end;

soil_albedo!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    sbulk = spac.soil_bulk;
    top_soil = spac.soils[1];

    @assert 1 <= sbulk.trait.color <=20;

    soil_albedo!(config, sbulk, top_soil, config.SOIL_ALBEDO);

    return nothing
);

soil_albedo!(config::SPACConfiguration{FT}, sbulk::SoilBulk{FT}, top_soil::SoilLayer{FT}, albedo::SoilAlbedoPrescribe) where {FT} = nothing;

soil_albedo!(config::SPACConfiguration{FT}, sbulk::SoilBulk{FT}, top_soil::SoilLayer{FT}, albedo::SoilAlbedoBroadbandCLM) where {FT} = (
    # use linear interpolation method or CLM method (with upper limit)
    delta = max(0, FT(0.11) - FT(0.4) * max(top_soil.trait.vc.Θ_SAT, top_soil.state.θ + top_soil.state.θ_ice));
    par::FT = min(SOIL_ALBEDOS[sbulk.trait.color,1], SOIL_ALBEDOS[sbulk.trait.color,3] + delta);
    nir::FT = min(SOIL_ALBEDOS[sbulk.trait.color,2], SOIL_ALBEDOS[sbulk.trait.color,4] + delta);

    soil_albedo!(config, sbulk, par, nir, false);

    return nothing
);

soil_albedo!(config::SPACConfiguration{FT}, sbulk::SoilBulk{FT}, top_soil::SoilLayer{FT}, albedo::SoilAlbedoHyperspectralCLM) where {FT} = (
    # use linear interpolation method or CLM method (with upper limit)
    delta = max(0, FT(0.11) - FT(0.4) * max(top_soil.trait.vc.Θ_SAT, top_soil.state.θ + top_soil.state.θ_ice));
    par::FT = min(SOIL_ALBEDOS[sbulk.trait.color,1], SOIL_ALBEDOS[sbulk.trait.color,3] + delta);
    nir::FT = min(SOIL_ALBEDOS[sbulk.trait.color,2], SOIL_ALBEDOS[sbulk.trait.color,4] + delta);

    soil_albedo!(config, sbulk, par, nir, true);

    return nothing
);

soil_albedo!(config::SPACConfiguration{FT}, sbulk::SoilBulk{FT}, top_soil::SoilLayer{FT}, albedo::SoilAlbedoBroadbandCLIMA) where {FT} = (
    # use linear interpolation method or CLM method (with upper limit)
    rwc = max(top_soil.trait.vc.Θ_SAT, top_soil.state.θ + top_soil.state.θ_ice) / top_soil.trait.vc.Θ_SAT;
    par::FT = SOIL_ALBEDOS[sbulk.trait.color,1] * (1 - rwc) + rwc * SOIL_ALBEDOS[sbulk.trait.color,3];
    nir::FT = SOIL_ALBEDOS[sbulk.trait.color,2] * (1 - rwc) + rwc * SOIL_ALBEDOS[sbulk.trait.color,4];

    soil_albedo!(config, sbulk, par, nir, false);

    return nothing
);

soil_albedo!(config::SPACConfiguration{FT}, sbulk::SoilBulk{FT}, top_soil::SoilLayer{FT}, albedo::SoilAlbedoHyperspectralCLIMA) where {FT} = (
    # use linear interpolation method or CLM method (with upper limit)
    rwc = max(top_soil.trait.vc.Θ_SAT, top_soil.state.θ + top_soil.state.θ_ice) / top_soil.trait.vc.Θ_SAT;
    par::FT = SOIL_ALBEDOS[sbulk.trait.color,1] * (1 - rwc) + rwc * SOIL_ALBEDOS[sbulk.trait.color,3];
    nir::FT = SOIL_ALBEDOS[sbulk.trait.color,2] * (1 - rwc) + rwc * SOIL_ALBEDOS[sbulk.trait.color,4];

    soil_albedo!(config, sbulk, par, nir, true);

    return nothing
);

soil_albedo!(config::SPACConfiguration{FT}, sbulk::SoilBulk{FT}, ρ_par::FT, ρ_nir::FT, hyperspectral::Bool) where {FT} = (
    (; SPECTRA) = config;

    # if not hyperspectral, use broadband
    if !hyperspectral
        sbulk.auxil.ρ_sw[SPECTRA.IΛ_PAR] .= ρ_par;
        sbulk.auxil.ρ_sw[SPECTRA.IΛ_NIR] .= ρ_nir;

        return nothing
    end;

    # if hyperspectral, use hyperspectral method
    # TODO: use a new soil moddel for this, do not use GSV which is not process-based
    ρ_sw = similar(sbulk.auxil.ρ_sw);
    ρ_sw[SPECTRA.IΛ_PAR] .= ρ_par;
    ρ_sw[SPECTRA.IΛ_NIR] .= ρ_nir;
    sbulk.auxil.weight .= pinv(SPECTRA.MAT_SOIL) * ρ_sw;

    # function to solve for weights
    @inline _fit(x::Vector{FT}) where {FT} = (
        mul!(ρ_sw, SPECTRA.MAT_SOIL, x);
        tmp_vec_nir = abs.(view(ρ_sw,SPECTRA.IΛ_NIR) .- ρ_nir);
        diff = ( mean( view(ρ_sw,SPECTRA.IΛ_PAR) ) - ρ_par ) ^ 2 + mean( tmp_vec_nir ) ^ 2;

        return -diff
    );

    # solve for weights
    ms = ReduceStepMethodND{FT}(x_mins = FT[-2,-2,-2,-2], x_maxs = FT[2,2,2,2], x_inis = sbulk.auxil.weight, Δ_inis = FT[0.1,0.1,0.1,0.1]);
    tol = SolutionToleranceND{FT}(FT[0.001,0.001,0.001,0.001], 50);
    sol = find_peak(_fit, ms, tol);
    sbulk.auxil.weight .= sol;

    # update vectors in soil
    mul!(sbulk.auxil.ρ_sw, SPECTRA.MAT_SOIL, sbulk.auxil.weight);

    # make sure the albedo is at least 0.01
    @. sbulk.auxil.ρ_sw = max.(sbulk.auxil.ρ_sw, FT(0.01));

    return nothing
);
