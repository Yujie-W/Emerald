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


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-02: generalize the function from CanopyLayers.dladgen to inclination_angles!
#     2022-Jun-02: add method for VerhoefLIDF algorithm
#     2023-May-22: add support to BetaLIDF
#     2023-Oct-14: if LAI <= 0, do nothing
#     2023-Oct-18: if LAI <= 0 && SAI <= 0, do nothing
#
#######################################################################################################################################################################################################
"""

    inclination_angles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Update the frequency of leaf inclination angles, given
- `config` SPAC configurations
- `spac` `BulkSPAC` type multiple layer SPAC

"""
function inclination_angles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    can_str = spac.canopy.structure;

    if can_str.state.lai <= 0 && can_str.state.sai <= 0
        return nothing
    end;

    (; Θ_INCL_BNDS) = config;

    # TODO: make p_incl an auxiliary variable
    for i in eachindex(can_str.state.p_incl)
        can_str.state.p_incl[i] = lidf_cdf(can_str.state.lidf, Θ_INCL_BNDS[i,2]) - lidf_cdf(can_str.state.lidf, Θ_INCL_BNDS[i,1]);
    end;

    return nothing
end;
