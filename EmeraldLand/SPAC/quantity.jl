#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Nov-18: add function to read beta from spac
#
#######################################################################################################################################################################################################
"""

    BETA(spac::BulkSPAC{FT}) where {FT}

Return the average beta factor for
- `spac` `BulkSPAC` SPAC

"""
function BETA end;

BETA(spac::BulkSPAC{FT}) where {FT} = (
    leaves = spac.plant.leaves;

    # compute the mean beta
    βs = 0;
    for leaf in leaves
        βs += read_β(leaf);
    end;

    return βs / length(leaves)
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-19: add function to compute canopy net primary productivity
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

    CNPP(spac::BulkSPAC{FT}) where {FT}

Return the canopy net primary productivity per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function CNPP end;

CNPP(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute GPP
    cnpp::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        cnpp += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_n_sunlit) +
                 (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_n_shaded) * canopy.structure.state.δlai[j];
    end;

    return cnpp
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Sep-07: add function to add up GPP for SPAC
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

    GPP(spac::BulkSPAC{FT}) where {FT}

Return the gross primary productivity per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function GPP end;

GPP(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute GPP
    gpp::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        gpp += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit) +
                (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded) * canopy.structure.state.δlai[j];
    end;

    return gpp
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-11: add function to compute PAR above canopy
#
#######################################################################################################################################################################################################
"""

    PAR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return the PAR above canopy per ground area, given
- `config` `SPACConfiguration` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function PAR end;

PAR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    (; SPECTRA) = config;
    rad_sw = spac.meteo.rad_sw;

    ppfd = photon.(SPECTRA.Λ_PAR, (rad_sw.e_dif + rad_sw.e_dir)[SPECTRA.IΛ_PAR]) .* 1000;

    return ppfd' * SPECTRA.ΔΛ_PAR
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-19: add function to compute canopy integrated PPAR
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

    PPAR(spac::BulkSPAC{FT}) where {FT}

Return the canopy integrated PPAR per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function PPAR end;

PPAR(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute GPP
    ppar::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        ppar += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.ppar_sunlit) +
                 (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.ppar_shaded) * canopy.structure.state.δlai[i];
    end;

    return ppar
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Sep-08: add function to add up transpiration rate for SPAC
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

    T_VEG(spac::BulkSPAC{FT}) where {FT}

Return the transpiration rate per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function T_VEG end;

T_VEG(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute transpiration rate
    tran::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        tran += flow_out(leaves[i]) * canopy.structure.state.δlai[j];
    end;

    return tran
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-09: add function to compute the weighted average of ϕ_d, ϕ_f, ϕ_n, and ϕ_p
#
#######################################################################################################################################################################################################

function ΦDFNP end;

ΦDFNP(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    sum_ϕda::FT = 0;
    sum_ϕfa::FT = 0;
    sum_ϕna::FT = 0;
    sum_ϕpa::FT = 0;
    sum_a::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        sum_ϕda += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit .* leaves[i].flux.auxil.ϕ_d_sunlit) +
                   (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded * leaves[i].flux.auxil.ϕ_d_shaded) * canopy.structure.state.δlai[j];
        sum_ϕfa += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit .* leaves[i].flux.auxil.ϕ_f_sunlit) +
                   (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded * leaves[i].flux.auxil.ϕ_f_shaded) * canopy.structure.state.δlai[j];
        sum_ϕna += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit .* leaves[i].flux.auxil.ϕ_n_sunlit) +
                   (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded * leaves[i].flux.auxil.ϕ_n_shaded) * canopy.structure.state.δlai[j];
        sum_ϕpa += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit .* leaves[i].flux.auxil.ϕ_p_sunlit) +
                   (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded * leaves[i].flux.auxil.ϕ_p_shaded) * canopy.structure.state.δlai[j];
        sum_a   += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit) +
                    (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded) * canopy.structure.state.δlai[j];
    end;

    return (sum_ϕda, sum_ϕfa, sum_ϕna, sum_ϕpa) ./ sum_a
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jun-13: add function to add up total ETR
#
#######################################################################################################################################################################################################
"""

    ΣETR(spac::BulkSPAC{FT}) where {FT}

Return the total ETR per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function ΣETR end;

ΣETR(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute GPP
    Σetr::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        Σetr += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.etr_sunlit) +
                 (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.etr_shaded) * canopy.structure.state.δlai[j];
    end;

    return Σetr
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jun-13: add function to add up total SIF photons (unit mol m⁻² s⁻¹)
#
#######################################################################################################################################################################################################
"""

ΣSIF(spac::BulkSPAC{FT}) where {FT}

Return the total SIF at chloroplast level (without any reabsorption) per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function ΣSIF end;

ΣSIF(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute SIF in photons unit
    Σsif::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        Σsif += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.ppar_sunlit .* leaves[i].flux.auxil.ϕ_f_sunlit) +
                 (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.ppar_shaded * leaves[i].flux.auxil.ϕ_f_shaded) * canopy.structure.state.δlai[j];
    end;

    return Σsif
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-11: add function to add up total SIF at chloroplast level (without any reabsorption)
#     2023-Sep-11: convert the unit to W m⁻² per ground area from nW m⁻²
#
#######################################################################################################################################################################################################
"""

    ΣSIF_CHL(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return the total SIF at chloroplast level (without any reabsorption) in W m⁻² per ground area, given
- `config` `SPACConfiguration` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function ΣSIF_CHL end;

ΣSIF_CHL(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    (; SPECTRA) = config;
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute SIF in energy unit before reabsorption within leaves (W m⁻²)
    Σsif::FT = 0;
    for i in eachindex(leaves)
        Σsif += view(canopy.sun_geometry.auxil.e_sif_chl,:,i)' * SPECTRA.ΔΛ_SIF / 1000;
    end;

    return Σsif
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-11: add function to add up total SIF in energy unit at leaf level after reabsorption
#     2023-Sep-11: convert the unit to W m⁻² per ground area from nW m⁻²
#
#######################################################################################################################################################################################################
"""

    ΣSIF_LEAF(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return the total SIF at leaf level after reabsorption in W m⁻² per ground area, given
- `config` `SPACConfiguration` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function ΣSIF_LEAF end;

ΣSIF_LEAF(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    (; SPECTRA) = config;
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute SIF in energy unit after reabsorption within leaves (W m⁻²)
    Σsif::FT = 0;
    for i in eachindex(leaves)
        Σsif += (canopy.sun_geometry.auxil.e_sifꜜ_layer[:,i] .+ canopy.sun_geometry.auxil.e_sifꜛ_layer[:,i])' * SPECTRA.ΔΛ_SIF / 1000;
    end;

    return Σsif
);
