# This file contains funtion to compute leaf inclination angle distribution

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
    θ = deg2rad(θ);
    y::FT = 0;
    x::FT = 2θ;
    n = 0;
    δx::FT = 1;
    while (abs(δx) >= max(eps(FT), 1e-8)) && (n < 50)
        y = A * sin(x) + B / 2 * sin(2x);
        δx = (y - x) / 2 + θ;
        x += δx;
        n += 1;
    end;

    return 2 * (y + θ) / FT(π)
);
