"""

    colimited_rate(a_1::FT, a_2::FT, colim::Union{MinimumColimit,QuadraticColimit{FT},SerialColimit,SquareColimit}) where {FT}

Return the minimum of two rates, given
- `a_1` Rate 1
- `a_2` Rate 2
- `colim` `MinimumColimit`, `QuadraticColimit`, `SerialColimit`, or `SquareColimit` type struct

"""
function colimited_rate end;

colimited_rate(a_1::FT, a_2::FT, ::MinimumColimit) where {FT} = min(a_1, a_2);

colimited_rate(a_1::FT, a_2::FT, colim::QuadraticColimit{FT}) where {FT} = lower_quadratic(colim.CURVATURE, -(a_1 + a_2), a_1 * a_2);

colimited_rate(a_1::FT, a_2::FT, ::SerialColimit) where {FT} = a_1 * a_2 / (a_1 + a_2);

colimited_rate(a_1::FT, a_2::FT, ::SquareColimit) where {FT} = a_1 * a_2 / sqrt(a_1 ^ 2 + a_2 ^ 2);


"""

    colimited_rate!(a_1::Union{FT,Vector{FT}}, a_2::Vector{FT}, a_i::Vector{FT}, colim::Union{MinimumColimit,QuadraticColimit{FT},SerialColimit,SquareColimit}) where {FT}

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
            ::MinimumColimit) where {FT} = (@. a_i = min(a_1, a_2); return nothing);

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
            ::SerialColimit) where {FT} = (@. a_i = a_1 * a_2 / (a_1 + a_2); return nothing);

colimited_rate!(
            a_1::Union{FT, Vector{FT}},
            a_2::Vector{FT},
            a_i::Vector{FT},
            ::SquareColimit) where {FT} = (@. a_i = a_1 * a_2 / sqrt(a_1 ^ 2 + a_2 ^ 2); return nothing);


"""

    colimit_photosynthesis!(psm::LeafPhotosystem{FT}; β::FT = FT(1)) where {FT}
    colimit_photosynthesis!(psm::CanopyLayerPhotosystem{FT}; β::FT = FT(1)) where {FT}

Colimit the photosynthesis by rubisco-, light-, and product-limited photosynthetic rates, given
- `psm` `CanopyLayerPhotosystem` or `LeafPhotosystem` type photosynthesis model
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd (default is 1)

"""
function colimit_photosynthesis! end;

colimit_photosynthesis!(config::SPACConfig{FT}, psm::CanopyLayerPhotosystem{FT}; β::FT = FT(1)) where {FT} = (
    colimited_rate!(psm.auxil.a_c, psm.auxil.a_j, psm.auxil.a_i, config.METHODS.COLIMIT_CJ);
    colimited_rate!(psm.auxil.a_p, psm.auxil.a_i, psm.auxil.a_g, config.METHODS.COLIMIT_IP);
    @. psm.auxil.a_n = psm.auxil.a_g - β .* psm.auxil.r_d;

    return nothing
);

colimit_photosynthesis!(config::SPACConfig{FT}, psm::LeafPhotosystem{FT}; β::FT = FT(1)) where {FT} = (
    a_i = colimited_rate(psm.auxil.a_c, psm.auxil.a_j, config.METHODS.COLIMIT_CJ);
    psm.auxil.a_g = colimited_rate(psm.auxil.a_p, a_i, config.METHODS.COLIMIT_IP);
    psm.auxil.a_n = psm.auxil.a_g - β * psm.auxil.r_d;

    return nothing
);
