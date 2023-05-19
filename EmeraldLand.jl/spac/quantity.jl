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
    (; DIM_LAYER, LEAVES) = spac;

    # compute the mean beta
    _βs = 0;
    for _leaves in LEAVES
        _βs += β_factor(_leaves.SM);
    end;

    return _βs / DIM_LAYER
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
    (; CANOPY, DIM_LAYER, LEAVES) = spac;

    # compute GPP
    _cnpp::FT = 0;
    for _i in 1:DIM_LAYER
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
    (; CANOPY, DIM_LAYER, LEAVES) = spac;

    # compute GPP
    _gpp::FT = 0;
    for _i in 1:DIM_LAYER
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
    (; CANOPY, DIM_LAYER, LEAVES) = spac;

    # compute GPP
    _ppar::FT = 0;
    for _i in 1:DIM_LAYER
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
    (; CANOPY, DIM_LAYER, LEAVES) = spac;

    # compute transpiration rate
    _tran::FT = 0;
    for _i in 1:DIM_LAYER
        _tran += flow_out(LEAVES[_i]) * CANOPY.δlai[_i];
    end;

    return _tran
);
