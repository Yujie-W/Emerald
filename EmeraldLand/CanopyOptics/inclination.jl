#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-02: generalize the function from CanopyLayers.dcum to lidf_cdf
#     2023-May-22: add method for BetaLIDF
#     2023-May-31: VerhoefLIDF judgement to A > 1 from A >= 1
#
#######################################################################################################################################################################################################
"""

    lidf_cdf(lidf::BetaLIDF{FT}, θ::FT) where {FT}
    lidf_cdf(lidf::VerhoefLIDF{FT}, θ::FT) where {FT}

Return the cumulative distribution frequency, given
- `lidf` `BetaLIDF` or `VerhoefLIDF` type algorithm
- `θ` Leaf inclination angle in `[°]`

"""
function lidf_cdf end;

lidf_cdf(lidf::BetaLIDF{FT}, θ::FT) where {FT} = (
    (; A, B) = lidf;

    return beta_inc(A, B, θ / 90, 1 - θ / 90)[1]
);

lidf_cdf(lidf::VerhoefLIDF{FT}, θ::FT) where {FT} = (
    (; A, B) = lidf;

    if A > 1
        return 1 - cosd(θ)
    end;

    # iterate to solve for CDF solution
    _θ = deg2rad(θ);
    _y::FT = 0;
    _x = 2 * _θ;
    _n = 0;
    _δx::FT = 1;
    while (abs(_δx) >= max(eps(FT), 1e-8)) && (_n < 50)
        _y = A * sin(_x) + B / 2 * sin(2*_x);
        _δx = (_y - _x) / 2 + _θ;
        _x += _δx;
        _n += 1;
    end;

    return 2 * (_y + _θ) / FT(π)
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-02: generalize the function from CanopyLayers.dladgen to inclination_angles!
#     2022-Jun-02: add method for VerhoefLIDF algorithm
#     2023-May-22: add sypport to BetaLIDF
#     2023-Jun-20: add config to parameter list
#
#######################################################################################################################################################################################################
"""

    inclination_angles!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}, lidf::Union{BetaLIDF{FT}, VerhoefLIDF{FT}}) where {FT}

Update the frequency of leaf inclination angles, given
- `config` SPAC configurations
- `can` `MultiLayerCanopy` type multiple layer canopy
- `lidf` `BetaLIDF` or `VerhoefLIDF` type algorithm

"""
function inclination_angles! end;

inclination_angles!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = inclination_angles!(config, spac.CANOPY, spac.CANOPY.LIDF);

inclination_angles!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}, lidf::Union{BetaLIDF{FT}, VerhoefLIDF{FT}}) where {FT} = (
    (; Θ_INCL_BNDS) = config;

    for i in eachindex(can.P_INCL)
        can.P_INCL[i] = lidf_cdf(lidf, Θ_INCL_BNDS[i,2]) - lidf_cdf(lidf, Θ_INCL_BNDS[i,1]);
    end;

    return nothing
);
