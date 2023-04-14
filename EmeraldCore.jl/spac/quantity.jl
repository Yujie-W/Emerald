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

BETA(spac::MultiLayerSPAC{FT,DIMS}) where {FT,DIMS} = (
    (; LEAVES) = spac;
    (; DIM_CANOPY) = DIMS;

    # compute the mean beta
    _βs = 0;
    for _leaves in LEAVES
        _βs += β_factor(_leaves.SM);
    end;

    return _βs / DIM_CANOPY
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-19: add function to compute canopy net primary productivity
#
#######################################################################################################################################################################################################
"""

    CNPP(spac::MultiLayerSPAC{FT}) where {FT}

Return the canopy net primary productivity per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function CNPP end

CNPP(spac::MultiLayerSPAC{FT,DIMS}) where {FT,DIMS} = (
    (; CANOPY, LEAVES) = spac;
    (; DIM_CANOPY) = DIMS;

    # compute GPP
    _cnpp::FT = 0;
    for _i in 1:DIM_CANOPY
        _cnpp += CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_net_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_net_shaded;
    end;
    _cnpp *= spac.CANOPY.lai / DIM_CANOPY;

    return _cnpp
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Sep-07: add function to add up GPP for SPAC
#
#######################################################################################################################################################################################################
"""

    GPP(spac::MultiLayerSPAC{FT}) where {FT}

Return the gross primary productivity per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function GPP end

GPP(spac::MultiLayerSPAC{FT,DIMS}) where {FT,DIMS} = (
    (; CANOPY, LEAVES) = spac;
    (; DIM_CANOPY) = DIMS;

    # compute GPP
    _gpp::FT = 0;
    for _i in 1:DIM_CANOPY
        _gpp += CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].a_gross_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].a_gross_shaded;
    end;
    _gpp *= spac.CANOPY.lai / DIM_CANOPY;

    return _gpp
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-19: add function to compute canopy integrated PPAR
#
#######################################################################################################################################################################################################
"""

    PPAR(spac::MultiLayerSPAC{FT}) where {FT}

Return the canopy integrated PPAR per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function PPAR end

PPAR(spac::MultiLayerSPAC{FT,DIMS}) where {FT,DIMS} = (
    (; CANOPY, LEAVES) = spac;
    (; DIM_CANOPY) = DIMS;

    # compute GPP
    _ppar::FT = 0;
    for _i in 1:DIM_CANOPY
        _ppar += CANOPY.OPTICS.p_sunlit[_i] * mean(LEAVES[_i].ppar_sunlit) + (1 - CANOPY.OPTICS.p_sunlit[_i]) * LEAVES[_i].ppar_shaded;
    end;
    _ppar *= spac.CANOPY.lai / DIM_CANOPY;

    return _ppar
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Sep-08: add function to add up transpiration rate for SPAC
#
#######################################################################################################################################################################################################
"""

    T_VEG(spac::MultiLayerSPAC{FT}) where {FT}

Return the transpiration rate per ground area, given
- `spac` `MultiLayerSPAC` SPAC

"""
function T_VEG end

T_VEG(spac::MultiLayerSPAC{FT,DIMS}) where {FT,DIMS} = (
    (; LEAVES) = spac;
    (; DIM_CANOPY) = DIMS;

    # compute transpiration rate
    _tran::FT = 0;
    for _i in 1:DIM_CANOPY
        _tran += flow_out(LEAVES[_i]);
    end;
    _tran *= spac.CANOPY.lai / DIM_CANOPY;

    return _tran
);
