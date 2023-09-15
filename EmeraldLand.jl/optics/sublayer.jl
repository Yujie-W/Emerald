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


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the absorption ratio from chlorophyll a and b
#
#######################################################################################################################################################################################################
"""

    sublayer_f_cab(bios::HyperLeafBioState{FT}, k_ant::FT, k_brown::FT, k_cab::FT, k_car_v::FT, k_car_z::FT, k_cbc::FT, k_H₂O::FT, k_lma::FT, k_pro::FT, lwc::FT) where {FT}

Return the fraction of chlorophyll a and b absorption in the sublayer, given
- `bios` leaf biophysical state variables
- `k_ant` anthocynanin absorption coefficient
- `k_brown` brown pigments absorption coefficient
- `k_cab` chlorophyll a and b absorption coefficient
- `k_car_v` violaxanthin carotenoid absorption coefficient
- `k_car_z` zeaxanthin carotenoid absorption coefficient
- `k_cbc` carbon-based constituents absorption coefficient
- `k_H₂O` water absorption coefficient
- `k_lma` dry mass absorption coefficient
- `k_pro` protein absorption coefficient
- `lwc` leaf water content

"""
function sublayer_f_cab(bios::HyperLeafBioState{FT}, k_ant::FT, k_brown::FT, k_cab::FT, k_car_v::FT, k_car_z::FT, k_cbc::FT, k_H₂O::FT, k_lma::FT, k_pro::FT, lwc::FT) where {FT}
    Σkx = k_ant * bios.ant +                            # anthocynanin absorption absorption
          k_brown * bios.brown +                        # brown pigments
          k_cab * bios.cab +                            # chlorophyll a + b absorption
          k_car_v * bios.car * (1 - bios.f_zeax) +      # violaxanthin carotenoid absorption
          k_car_z * bios.car * bios.f_zeax +            # zeaxanthin carotenoid absorption
          k_cbc * bios.cbc +                            # carbon-based constituents absorption
          k_H₂O * (lwc * M_H₂O(FT) / ρ_H₂O(FT) * 100) + # water absorption
          k_lma * (bios.lma - bios.cbc - bios.pro) +    # dry mass absorption (if some remained)
          k_pro * bios.pro;                             # protein absorption

    return k_cab * bios.cab / Σkx
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the absorption ratio from carotenoid
#
#######################################################################################################################################################################################################
"""

    sublayer_f_car(bios::HyperLeafBioState{FT}, k_ant::FT, k_brown::FT, k_cab::FT, k_car_v::FT, k_car_z::FT, k_cbc::FT, k_H₂O::FT, k_lma::FT, k_pro::FT, lwc::FT) where {FT}

Return the fraction of chlorophyll a and b absorption in the sublayer, given
- `bios` leaf biophysical state variables
- `k_ant` anthocynanin absorption coefficient
- `k_brown` brown pigments absorption coefficient
- `k_cab` chlorophyll a and b absorption coefficient
- `k_car_v` violaxanthin carotenoid absorption coefficient
- `k_car_z` zeaxanthin carotenoid absorption coefficient
- `k_cbc` carbon-based constituents absorption coefficient
- `k_H₂O` water absorption coefficient
- `k_lma` dry mass absorption coefficient
- `k_pro` protein absorption coefficient
- `lwc` leaf water content

"""
function sublayer_f_car(bios::HyperLeafBioState{FT}, k_ant::FT, k_brown::FT, k_cab::FT, k_car_v::FT, k_car_z::FT, k_cbc::FT, k_H₂O::FT, k_lma::FT, k_pro::FT, lwc::FT) where {FT}
    Σkx = k_ant * bios.ant +
          k_brown * bios.brown +
          k_cab * bios.cab +
          k_car_v * bios.car * (1 - bios.f_zeax) +
          k_car_z * bios.car * bios.f_zeax +
          k_cbc * bios.cbc +
          k_H₂O * (lwc * M_H₂O(FT) / ρ_H₂O(FT) * 100) +
          k_lma * (bios.lma - bios.cbc - bios.pro) +
          k_pro * bios.pro;

    return (k_car_v * bios.car * (1 - bios.f_zeax) + k_car_z * bios.car * bios.f_zeax) / Σkx
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the sublayer transmission
#
#######################################################################################################################################################################################################
"""

    sublayer_τ(bios::HyperLeafBioState{FT}, k_ant::FT, k_brown::FT, k_cab::FT, k_car_v::FT, k_car_z::FT, k_cbc::FT, k_H₂O::FT, k_lma::FT, k_pro::FT, lwc::FT, x::FT, N::Int) where {FT}

Return the fraction of chlorophyll a and b absorption in the sublayer, given
- `bios` leaf biophysical state variables
- `k_ant` anthocynanin absorption coefficient
- `k_brown` brown pigments absorption coefficient
- `k_cab` chlorophyll a and b absorption coefficient
- `k_car_v` violaxanthin carotenoid absorption coefficient
- `k_car_z` zeaxanthin carotenoid absorption coefficient
- `k_cbc` carbon-based constituents absorption coefficient
- `k_H₂O` water absorption coefficient
- `k_lma` dry mass absorption coefficient
- `k_pro` protein absorption coefficient
- `lwc` leaf water content
- `x` proportion of the layer of the whole leaf
- `N` number of sublayers of the whole layer

"""
function sublayer_τ(bios::HyperLeafBioState{FT}, k_ant::FT, k_brown::FT, k_cab::FT, k_car_v::FT, k_car_z::FT, k_cbc::FT, k_H₂O::FT, k_lma::FT, k_pro::FT, lwc::FT, x::FT, N::Int) where {FT}
    Σkx = k_ant * bios.ant +
          k_brown * bios.brown +
          k_cab * bios.cab +
          k_car_v * bios.car * (1 - bios.f_zeax) +
          k_car_z * bios.car * bios.f_zeax +
          k_cbc * bios.cbc +
          k_H₂O * (lwc * M_H₂O(FT) / ρ_H₂O(FT) * 100) +
          k_lma * (bios.lma - bios.cbc - bios.pro) +
          k_pro * bios.pro;
    Σkx *= x;
    τ_layer = (1 - Σkx) * exp(-Σkx) + Σkx^2 * expint(Σkx + eps(FT));

    return τ_layer ^ (FT(1) / N)
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to update the sublayer absorption and transmittance within HyperLeafBio
#
#######################################################################################################################################################################################################
"""

    leaf_sublayer_absorptions!(config::SPACConfiguration{FT}, bio::HyperLeafBio{FT}, lwc::FT; N::Int = 10) where {FT}
    leaf_sublayer_absorptions!(lha::HyperspectralAbsorption{FT}, bio::HyperLeafBio{FT}, lwc::FT; N::Int = 10) where {FT}

Update the sublayer absorption and transmittance within `bio`, given
- `config` SPAC configuration
- `bio` HyperLeafBio struct
- `lwc` leaf water content
- `N` number of sublayers of the whole layer
- `lha` HyperspectralAbsorption struct

"""
function leaf_sublayer_absorptions! end;

leaf_sublayer_absorptions!(config::SPACConfiguration{FT}, bio::HyperLeafBio{FT}, lwc::FT; N::Int = 10) where {FT} = leaf_sublayer_absorptions!(config.lha, bio, lwc; N = N);

leaf_sublayer_absorptions!(lha::HyperspectralAbsorption{FT}, bio::HyperLeafBio{FT}, lwc::FT; N::Int = 10) where {FT} = (
    (; K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO) = lha;

    x = 1 / bio.state.meso_n;
    bio.auxil.f_cab .= sublayer_f_cab.((bio.state,), K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, lwc);
    bio.auxil.f_car .= sublayer_f_car.((bio.state,), K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, lwc);

    bio.auxil.τ_sub_1 .= sublayer_τ.((bio.state,), K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, lwc, x, N);
    bio.auxil.τ_sub_2 .= sublayer_τ.((bio.state,), K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, lwc, 1-x, N);

    return nothing
);
