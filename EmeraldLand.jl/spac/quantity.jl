#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Nov-18: add function to read beta from spac
#
#######################################################################################################################################################################################################
"""

    BETA(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Return the average beta factor for
- `spac` `MultiLayerSPAC` SPAC

"""
function BETA end

BETA(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; LEAVES) = spac;

    # compute the mean beta
    _βs = 0;
    for _leaves in LEAVES
        _βs += β_factor(_leaves.SM);
    end;

    return _βs / length(LEAVES)
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

    CNPP(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Return the canopy net primary productivity per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function CNPP end

CNPP(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; CANOPY, LEAVES) = spac;

    # compute GPP
    _cnpp::FT = 0;
    for _i in eachindex(LEAVES)
        _cnpp += (CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_net_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_net_shaded) * CANOPY.δlai[_i];
    end;

    return _cnpp
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

    GPP(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Return the gross primary productivity per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function GPP end

GPP(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; CANOPY, LEAVES) = spac;

    # compute GPP
    _gpp::FT = 0;
    for _i in eachindex(LEAVES)
        _gpp += (CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_gross_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_gross_shaded) * CANOPY.δlai[_i];
    end;

    return _gpp
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

    PPAR(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Return the canopy integrated PPAR per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function PPAR end

PPAR(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; CANOPY, LEAVES) = spac;

    # compute GPP
    _ppar::FT = 0;
    for _i in eachindex(LEAVES)
        _ppar += (CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].ppar_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].ppar_shaded) * CANOPY.δlai[_i];
    end;

    return _ppar
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

    T_VEG(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Return the transpiration rate per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function T_VEG end

T_VEG(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; CANOPY, LEAVES) = spac;

    # compute transpiration rate
    _tran::FT = 0;
    for _i in eachindex(LEAVES)
        _tran += flow_out(LEAVES[_i]) * CANOPY.δlai[_i];
    end;

    return _tran
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-09: add function to compute the weighted average of ϕ_d, ϕ_f, ϕ_n, and ϕ_p
#
#######################################################################################################################################################################################################

function ΦDFNP end

ΦDFNP(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; CANOPY, LEAVES) = spac;

    _sum_ϕda::FT = 0;
    _sum_ϕfa::FT = 0;
    _sum_ϕna::FT = 0;
    _sum_ϕpa::FT = 0;
    _sum_a::FT = 0;
    for _i in eachindex(LEAVES)
        _sum_ϕda += (CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_gross_sunlit .* LEAVES[_i].ϕ_d_sunlit) +
                    (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_gross_shaded * LEAVES[_i].ϕ_d_shaded) * CANOPY.δlai[_i];
        _sum_ϕfa += (CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_gross_sunlit .* LEAVES[_i].ϕ_f_sunlit) +
                    (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_gross_shaded * LEAVES[_i].ϕ_f_shaded) * CANOPY.δlai[_i];
        _sum_ϕna += (CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_gross_sunlit .* LEAVES[_i].ϕ_n_sunlit) +
                    (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_gross_shaded * LEAVES[_i].ϕ_n_shaded) * CANOPY.δlai[_i];
        _sum_ϕpa += (CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_gross_sunlit .* LEAVES[_i].ϕ_p_sunlit) +
                    (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_gross_shaded * LEAVES[_i].ϕ_p_shaded) * CANOPY.δlai[_i];
        _sum_a += (CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_gross_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_gross_shaded) * CANOPY.δlai[_i];
    end;

    return (_sum_ϕda, _sum_ϕfa, _sum_ϕna, _sum_ϕpa) ./ _sum_a
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jun-13: add function to add up total ETR
#
#######################################################################################################################################################################################################
"""

    ΣETR(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Return the total ETR per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function ΣETR end

ΣETR(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; CANOPY, LEAVES) = spac;

    # compute GPP
    _Σetr::FT = 0;
    for _i in eachindex(LEAVES)
        _Σetr += (CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].etr_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].etr_shaded) * CANOPY.δlai[_i];
    end;

    return _Σetr
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jun-13: add function to add up total SIF photons (unit mol m⁻² s⁻¹)
#
#######################################################################################################################################################################################################
"""

ΣSIF(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Return the total SIF at chloroplast level (without any reabsorption) per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function ΣSIF end

ΣSIF(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; CANOPY, LEAVES) = spac;

    # compute SIF in photons unit
    _Σsif::FT = 0;
    for _i in eachindex(LEAVES)
        _Σsif += (CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].ppar_sunlit .* LEAVES[_i].ϕ_f_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].ppar_shaded * LEAVES[_i].ϕ_f_shaded) *
                 CANOPY.δlai[_i];
    end;

    return _Σsif
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

    ΣSIF_CHL(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Return the total SIF at chloroplast level (without any reabsorption) in W m⁻² per ground area, given
- `config` `SPACConfiguration` SPAC configuration
- `spac` `MultiLayerSPAC` SPAC

"""
function ΣSIF_CHL end

ΣSIF_CHL(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; WLSET) = config;
    (; CANOPY, LEAVES) = spac;

    # compute SIF in energy unit before reabsorption within leaves (W m⁻²)
    _Σsif::FT = 0;
    for _i in eachindex(LEAVES)
        _Σsif += (CANOPY.RADIATION.s_layer_down_chl[:,_i] .+ CANOPY.RADIATION.s_layer_up_chl[:,_i])' * WLSET.ΔΛ_SIF / 1000;
    end;

    return _Σsif
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

    ΣSIF_LEAF(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Return the total SIF at leaf level after reabsorption in W m⁻² per ground area, given
- `config` `SPACConfiguration` SPAC configuration
- `spac` `MultiLayerSPAC` SPAC

"""
function ΣSIF_LEAF end

ΣSIF_LEAF(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; WLSET) = config;
    (; CANOPY, LEAVES) = spac;

    # compute SIF in energy unit after reabsorption within leaves (W m⁻²)
    _Σsif::FT = 0;
    for _i in eachindex(LEAVES)
        _Σsif += (CANOPY.RADIATION.s_layer_down[:,_i] .+ CANOPY.RADIATION.s_layer_up[:,_i])' * WLSET.ΔΛ_SIF / 1000;
    end;

    return _Σsif
);
