#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Nov-18: add function to read beta from spac
#
#######################################################################################################################################################################################################
"""

    BETA(spac::MultiLayerSPAC{FT}) where {FT}

Return the average beta factor for
- `spac` `MultiLayerSPAC` SPAC

"""
function BETA end

BETA(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; LEAVES) = spac;

    # compute the mean beta
    βs = 0;
    for _leaves in LEAVES
        βs += read_β(_leaves.flux.state.stomatal_model);
    end;

    return βs / length(LEAVES)
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

    CNPP(spac::MultiLayerSPAC{FT}) where {FT}

Return the canopy net primary productivity per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function CNPP end

CNPP(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; CANOPY, LEAVES) = spac;

    # compute GPP
    cnpp::FT = 0;
    N = length(LEAVES);
    for i in eachindex(LEAVES)
        j = N - i + 1;
        cnpp += (CANOPY.sun_geometry.auxil.p_sunlit[j] * mean(LEAVES[i].flux.auxil.a_n_sunlit) + (1 - CANOPY.sun_geometry.auxil.p_sunlit[j]) * LEAVES[i].flux.auxil.a_n_shaded) * CANOPY.δlai[j];
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

    GPP(spac::MultiLayerSPAC{FT}) where {FT}

Return the gross primary productivity per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function GPP end

GPP(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; CANOPY, LEAVES) = spac;

    # compute GPP
    gpp::FT = 0;
    N = length(LEAVES);
    for i in eachindex(LEAVES)
        j = N - i + 1;
        gpp += (CANOPY.sun_geometry.auxil.p_sunlit[j] * mean(LEAVES[i].flux.auxil.a_g_sunlit) + (1 - CANOPY.sun_geometry.auxil.p_sunlit[j]) * LEAVES[i].flux.auxil.a_g_shaded) * CANOPY.δlai[j];
    end;

    return gpp
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

    PPAR(spac::MultiLayerSPAC{FT}) where {FT}

Return the canopy integrated PPAR per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function PPAR end

PPAR(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; CANOPY, LEAVES) = spac;

    # compute GPP
    ppar::FT = 0;
    N = length(LEAVES);
    for i in eachindex(LEAVES)
        j = N - i + 1;
        ppar += (CANOPY.sun_geometry.auxil.p_sunlit[j] * mean(LEAVES[i].flux.auxil.ppar_sunlit) + (1 - CANOPY.sun_geometry.auxil.p_sunlit[j]) * LEAVES[i].flux.auxil.ppar_shaded) * CANOPY.δlai[i];
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

    T_VEG(spac::MultiLayerSPAC{FT}) where {FT}

Return the transpiration rate per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function T_VEG end

T_VEG(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; CANOPY, LEAVES) = spac;

    # compute transpiration rate
    tran::FT = 0;
    N = length(LEAVES);
    for i in eachindex(LEAVES)
        j = N - i + 1;
        tran += flow_out(LEAVES[i]) * CANOPY.δlai[j];
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

function ΦDFNP end

ΦDFNP(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; CANOPY, LEAVES) = spac;

    sum_ϕda::FT = 0;
    sum_ϕfa::FT = 0;
    sum_ϕna::FT = 0;
    sum_ϕpa::FT = 0;
    sum_a::FT = 0;
    N = length(LEAVES);
    for i in eachindex(LEAVES)
        j = N - i + 1;
        sum_ϕda += (CANOPY.sun_geometry.auxil.p_sunlit[j] * mean(LEAVES[i].flux.auxil.a_g_sunlit .* LEAVES[i].flux.auxil.ϕ_d_sunlit) +
                   (1 - CANOPY.sun_geometry.auxil.p_sunlit[j]) * LEAVES[i].flux.auxil.a_g_shaded * LEAVES[i].flux.auxil.ϕ_d_shaded) * CANOPY.δlai[j];
        sum_ϕfa += (CANOPY.sun_geometry.auxil.p_sunlit[j] * mean(LEAVES[i].flux.auxil.a_g_sunlit .* LEAVES[i].flux.auxil.ϕ_f_sunlit) +
                   (1 - CANOPY.sun_geometry.auxil.p_sunlit[j]) * LEAVES[i].flux.auxil.a_g_shaded * LEAVES[i].flux.auxil.ϕ_f_shaded) * CANOPY.δlai[j];
        sum_ϕna += (CANOPY.sun_geometry.auxil.p_sunlit[j] * mean(LEAVES[i].flux.auxil.a_g_sunlit .* LEAVES[i].flux.auxil.ϕ_n_sunlit) +
                   (1 - CANOPY.sun_geometry.auxil.p_sunlit[j]) * LEAVES[i].flux.auxil.a_g_shaded * LEAVES[i].flux.auxil.ϕ_n_shaded) * CANOPY.δlai[j];
        sum_ϕpa += (CANOPY.sun_geometry.auxil.p_sunlit[j] * mean(LEAVES[i].flux.auxil.a_g_sunlit .* LEAVES[i].flux.auxil.ϕ_p_sunlit) +
                   (1 - CANOPY.sun_geometry.auxil.p_sunlit[j]) * LEAVES[i].flux.auxil.a_g_shaded * LEAVES[i].flux.auxil.ϕ_p_shaded) * CANOPY.δlai[j];
        sum_a   += (CANOPY.sun_geometry.auxil.p_sunlit[j] * mean(LEAVES[i].flux.auxil.a_g_sunlit) + (1 - CANOPY.sun_geometry.auxil.p_sunlit[j]) * LEAVES[i].flux.auxil.a_g_shaded) * CANOPY.δlai[j];
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

    ΣETR(spac::MultiLayerSPAC{FT}) where {FT}

Return the total ETR per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function ΣETR end

ΣETR(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; CANOPY, LEAVES) = spac;

    # compute GPP
    Σetr::FT = 0;
    N = length(LEAVES);
    for i in eachindex(LEAVES)
        j = N - i + 1;
        Σetr += (CANOPY.sun_geometry.auxil.p_sunlit[j] * mean(LEAVES[i].flux.auxil.etr_sunlit) + (1 - CANOPY.sun_geometry.auxil.p_sunlit[j]) * LEAVES[i].flux.auxil.etr_shaded) * CANOPY.δlai[j];
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

ΣSIF(spac::MultiLayerSPAC{FT}) where {FT}

Return the total SIF at chloroplast level (without any reabsorption) per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function ΣSIF end

ΣSIF(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; CANOPY, LEAVES) = spac;

    # compute SIF in photons unit
    Σsif::FT = 0;
    N = length(LEAVES);
    for i in eachindex(LEAVES)
        j = N - i + 1;
        Σsif += (CANOPY.sun_geometry.auxil.p_sunlit[j] * mean(LEAVES[i].flux.auxil.ppar_sunlit .* LEAVES[i].flux.auxil.ϕ_f_sunlit) +
                 (1 - CANOPY.sun_geometry.auxil.p_sunlit[j]) * LEAVES[i].flux.auxil.ppar_shaded * LEAVES[i].flux.auxil.ϕ_f_shaded) * CANOPY.δlai[j];
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

    ΣSIF_CHL(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Return the total SIF at chloroplast level (without any reabsorption) in W m⁻² per ground area, given
- `config` `SPACConfiguration` SPAC configuration
- `spac` `MultiLayerSPAC` SPAC

"""
function ΣSIF_CHL end

ΣSIF_CHL(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; SPECTRA) = config;
    (; CANOPY, LEAVES) = spac;

    # compute SIF in energy unit before reabsorption within leaves (W m⁻²)
    Σsif::FT = 0;
    for i in eachindex(LEAVES)
        Σsif += (CANOPY.RADIATION.s_layer_down_chl[:,i] .+ CANOPY.RADIATION.s_layer_up_chl[:,i])' * SPECTRA.ΔΛ_SIF / 1000;
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

    ΣSIF_LEAF(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Return the total SIF at leaf level after reabsorption in W m⁻² per ground area, given
- `config` `SPACConfiguration` SPAC configuration
- `spac` `MultiLayerSPAC` SPAC

"""
function ΣSIF_LEAF end

ΣSIF_LEAF(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; SPECTRA) = config;
    (; CANOPY, LEAVES) = spac;

    # compute SIF in energy unit after reabsorption within leaves (W m⁻²)
    Σsif::FT = 0;
    for i in eachindex(LEAVES)
        Σsif += (CANOPY.RADIATION.s_layer_down[:,i] .+ CANOPY.RADIATION.s_layer_up[:,i])' * SPECTRA.ΔΛ_SIF / 1000;
    end;

    return Σsif
);
