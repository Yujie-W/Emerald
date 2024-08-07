# This file contains function to calculate colimitation of photosynthesis

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add colimit function back
#     2022-Feb-07: add C3Cytochrome to method (colimit j_p680 and j_p700 as well)
#     2022-Feb-11: add colimited_rate for general purpose in ETR as well as a_gross
#     2022-Mar-01: add colimit method for serial colimitation
#     2022-Mar-01: add colimit method for square colimitation
#
#######################################################################################################################################################################################################
"""

    colimited_rate(a_1::FT, a_2::FT, colim::MinimumColimit{FT}) where {FT}
    colimited_rate(a_1::FT, a_2::FT, colim::QuadraticColimit{FT}) where {FT}
    colimited_rate(a_1::FT, a_2::FT, colim::SerialColimit{FT}) where {FT}
    colimited_rate(a_1::FT, a_2::FT, colim::SquareColimit{FT}) where {FT}

Return the minimum of two rates, given
- `a_1` Rate 1
- `a_2` Rate 2
- `colim` `MinimumColimit`, `QuadraticColimit`, or `SerialColimit` type struct

"""
function colimited_rate end;

colimited_rate(a_1::FT, a_2::FT, colim::MinimumColimit{FT}) where {FT} = min(a_1, a_2);

colimited_rate(a_1::FT, a_2::FT, colim::QuadraticColimit{FT}) where {FT} = lower_quadratic(colim.CURVATURE, -(a_1 + a_2), a_1 * a_2);

colimited_rate(a_1::FT, a_2::FT, colim::SerialColimit{FT}) where {FT} = a_1 * a_2 / (a_1 + a_2);

colimited_rate(a_1::FT, a_2::FT, colim::SquareColimit{FT}) where {FT} = a_1 * a_2 / sqrt(a_1 ^ 2 + a_2 ^ 2);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Aug-01: add function colimited_rate!
#
#######################################################################################################################################################################################################
"""

    colimited_rate!(a_1::Union{FT, Vector{FT}}, a_2::Vector{FT}, a_i::Vector{FT}, colim) where {FT}

Colimit the rates, given
- `a_1` Rate 1
- `a_2` Rate 2
- `a_i` Intermediate rate (overwritten)
- `colim` `MinimumColimit`, `QuadraticColimit`, `SerialColimit`, or `SquareColimit` type struct

"""
function colimited_rate! end;

colimited_rate!(
            a_1::Union{FT, Vector{FT}},
            a_2::Vector{FT},
            a_i::Vector{FT},
            colim::MinimumColimit{FT}) where {FT} = (@. a_i = min(a_1, a_2); return nothing);

# a_i .= lower_quadratic.(colim.CURVATURE, -a_1 .- a_2, a_1 .* a_2);
colimited_rate!(
            a_1::Union{FT, Vector{FT}},
            a_2::Vector{FT},
            a_i::Vector{FT},
            colim::QuadraticColimit{FT}) where {FT} = (@. a_i = lower_quadratic(colim.CURVATURE, -a_1 - a_2, a_1 * a_2); return nothing);

colimited_rate!(
            a_1::Union{FT, Vector{FT}},
            a_2::Vector{FT},
            a_i::Vector{FT},
            colim::SerialColimit{FT}) where {FT} = (@. a_i = a_1 * a_2 / (a_1 + a_2); return nothing);

colimited_rate!(
            a_1::Union{FT, Vector{FT}},
            a_2::Vector{FT},
            a_i::Vector{FT},
            colim::SquareColimit{FT}) where {FT} = (@. a_i = a_1 * a_2 / sqrt(a_1 ^ 2 + a_2 ^ 2); return nothing);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add colimit function back
#     2022-Jan-24: use colimit from psm to abstractize the MinimumColimit and QuadraticColimit methods
#     2022-Feb-07: add C3Cyto support
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#
#######################################################################################################################################################################################################
"""

    colimit_photosynthesis!(psm::LeafPhotosystem{FT}; β::FT = FT(1)) where {FT}

Colimit the photosynthesis by rubisco-, light-, and product-limited photosynthetic rates, given
- `psm` `LeafPhotosystem` type photosynthesis model
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd (default is 1)

"""
function colimit_photosynthesis! end;

colimit_photosynthesis!(psm::CanopyLayerPhotosystem{FT}; β::FT = FT(1)) where {FT} = (
    colimited_rate!(psm.auxil.a_c, psm.auxil.a_j, psm.auxil.a_i, psm.trait.COLIMIT_CJ);
    colimited_rate!(psm.auxil.a_p, psm.auxil.a_i, psm.auxil.a_g, psm.trait.COLIMIT_IP);
    @. psm.auxil.a_n = psm.auxil.a_g - β .* psm.auxil.r_d;

    return nothing
);

colimit_photosynthesis!(psm::LeafPhotosystem{FT}; β::FT = FT(1)) where {FT} = (
    a_i = colimited_rate(psm.auxil.a_c, psm.auxil.a_j, psm.trait.COLIMIT_CJ);
    psm.auxil.a_g = colimited_rate(psm.auxil.a_p, a_i, psm.trait.COLIMIT_IP);
    psm.auxil.a_n = psm.auxil.a_g - β * psm.auxil.r_d;

    return nothing
);
