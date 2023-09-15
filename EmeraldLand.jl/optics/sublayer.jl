#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-14: add function to compute the transmittance of a sublayer of a leaf with pigments
#     2023-Sep-14: compute the total τ for N layers and then downscale it to τ^(1/N)
#     2023-Sep-14: compute the fraction of cab and car absorption in the sublayer
#
#######################################################################################################################################################################################################
"""

    sublayer_τ(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, x::FT, N::Int) where {FT}

Return the transmittance of a sublayer of a leaf with pigments, given
- `lha` leaf hyperspectral absorption coefficients
- `bio` leaf hyperspectral biophysics
- `lwc` leaf water content
- `x` proportion of the layer of the whole leaf
- `N` number of sublayers of the whole layer

"""
function sublayer_τ(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, x::FT, N::Int) where {FT}
    (; K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO) = lha;

    # define the vectors
    _sum_ki_xi = similar(K_CAB);
    _k_cab = similar(K_CAB);
    _k_car = similar(K_CAB);
    _f_cab = similar(K_CAB);
    _f_car = similar(K_CAB);
    _τ_all = similar(K_CAB);

    # compute the sum of absorption coefficients in one of the N sublayers of 1 layer (x of the total leaf thickness)
    @. _sum_ki_xi = K_CAB   * bio.cab +                                 # chlorophyll absorption
                    K_CAR_V * bio.car * (1 - bio.f_zeax) +              # violaxanthin carotenoid absorption
                    K_CAR_Z * bio.car * bio.f_zeax +                    # zeaxanthin carotenoid absorption
                    K_ANT   * bio.ant +                                 # anthocynanin absorption absorption
                    K_BROWN * bio.brown +                               # brown pigments
                    K_H₂O   * (lwc * M_H₂O(FT) / ρ_H₂O(FT) * 100) +     # water absorption
                    K_CBC   * bio.cbc +                                 # carbon-based constituents absorption
                    K_PRO   * bio.pro +                                 # protein absorption
                    K_LMA   * (bio.lma - bio.cbc - bio.pro);            # dry mass absorption (if some remained)
    @. _k_cab = K_CAB * bio.cab;
    @. _k_car = K_CAR_V * bio.car * (1 - bio.f_zeax) + K_CAR_Z * bio.car * bio.f_zeax;

    # compute the absorption ratios from cab and car
    @. _f_cab = _k_cab / _sum_ki_xi;
    @. _f_car = _k_car / _sum_ki_xi;

    #_sum_ki_xi .*= x / N;
    _sum_ki_xi .*= x;

    # this is the case when light penetrate with an angle (integrated over all angles, without accounting for F_CELL)
    @. _τ_all = (1 - _sum_ki_xi) * exp(-_sum_ki_xi) + _sum_ki_xi^2 * expint(_sum_ki_xi + eps(FT));

    return _τ_all .^ (FT(1) / N), _f_cab, _f_car
end;
