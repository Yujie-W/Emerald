#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-27: add function to initialize SPAC
#     2022-Jun-27: add leaf area controller to make sure soil and leaf areas are consistent with leaf area index
#     2023-Mar-27: initialize soil and leaf e as well (because T, SWC may be changed)
#
#######################################################################################################################################################################################################
"""

    initialize!(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Initialize the SPAC, given
- `spac` `MultiLayerSPAC` SPAC

"""
function initialize! end

initialize!(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; CANOPY, DIM_LAYER, LEAVES, SOIL) = spac;

    # make sure soil energy is correctly scaled with temperature and soil water content
    for _slayer in SOIL.LAYERS
        _slayer.e = (_slayer.CP * _slayer.ρ + _slayer.θ * CP_L() * ρ_H₂O()) * _slayer.t;
    end;

    # make sure leaf area index setup and energy are correct
    for _clayer in LEAVES
        _clayer.HS.AREA = SOIL.AREA * CANOPY.lai / DIM_LAYER;
        _clayer.e = (_clayer.CP * _clayer.BIO.lma * 10 + _clayer.HS.v_storage * CP_L_MOL(FT)) * _clayer.t;
    end;

    # initialize leaf level spectra
    leaf_spectra!(spac);

    # initialize stomatal conductance
    stomatal_conductance!(spac, FT(0));

    return nothing
);
