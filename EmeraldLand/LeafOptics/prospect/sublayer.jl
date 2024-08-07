#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the absorption ratio from chlorophyll a and b
#
#######################################################################################################################################################################################################
"""

    sublayer_f_cab(biot::LeafBioTrait{FT}, bios::LeafBioState{FT}, k_ant::FT, k_brown::FT, k_cab::FT, k_car_v::FT, k_car_z::FT, k_cbc::FT, k_H₂O::FT, k_lma::FT, k_pro::FT, lwc::FT) where {FT}

Return the fraction of chlorophyll a and b absorption in the sublayer, given
- `biot` leaf biophysical trait variables
- `bios` leaf biophysical state variables
- `k_ant` anthocyanin absorption coefficient
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
function sublayer_f_cab(biot::LeafBioTrait{FT}, bios::LeafBioState{FT}, k_ant::FT, k_brown::FT, k_cab::FT, k_car_v::FT, k_car_z::FT, k_cbc::FT, k_H₂O::FT, k_lma::FT, k_pro::FT, lwc::FT) where {FT}
    Σkx = k_ant * biot.ant +                            # anthocyanin absorption absorption
          k_brown * biot.brown +                        # brown pigments
          k_cab * biot.cab +                            # chlorophyll a + b absorption
          k_car_v * biot.car * (1 - bios.f_zeax) +      # violaxanthin carotenoid absorption
          k_car_z * biot.car * bios.f_zeax +            # zeaxanthin carotenoid absorption
          k_cbc * biot.cbc +                            # carbon-based constituents absorption
          k_H₂O * (lwc * M_H₂O(FT) / ρ_H₂O(FT) * 100) + # water absorption
          k_lma * (biot.lma - biot.cbc - biot.pro) +    # dry mass absorption (if some remained)
          k_pro * biot.pro;                             # protein absorption

    return k_cab * biot.cab / Σkx
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the absorption ratio from carotenoid
#
#######################################################################################################################################################################################################
"""

    sublayer_f_car(biot::LeafBioTrait{FT}, bios::LeafBioState{FT}, k_ant::FT, k_brown::FT, k_cab::FT, k_car_v::FT, k_car_z::FT, k_cbc::FT, k_H₂O::FT, k_lma::FT, k_pro::FT, lwc::FT) where {FT}

Return the fraction of chlorophyll a and b absorption in the sublayer, given
- `bios` leaf biophysical state variables
- `bios` leaf biophysical state variables
- `k_ant` anthocyanin absorption coefficient
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
function sublayer_f_car(biot::LeafBioTrait{FT}, bios::LeafBioState{FT}, k_ant::FT, k_brown::FT, k_cab::FT, k_car_v::FT, k_car_z::FT, k_cbc::FT, k_H₂O::FT, k_lma::FT, k_pro::FT, lwc::FT) where {FT}
    Σkx = k_ant * biot.ant +
          k_brown * biot.brown +
          k_cab * biot.cab +
          k_car_v * biot.car * (1 - bios.f_zeax) +
          k_car_z * biot.car * bios.f_zeax +
          k_cbc * biot.cbc +
          k_H₂O * (lwc * M_H₂O(FT) / ρ_H₂O(FT) * 100) +
          k_lma * (biot.lma - biot.cbc - biot.pro) +
          k_pro * biot.pro;

    return (k_car_v * biot.car * (1 - bios.f_zeax) + k_car_z * biot.car * bios.f_zeax) / Σkx
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the sublayer transmission
#
#######################################################################################################################################################################################################
"""

    sublayer_τ(biot::LeafBioTrait{FT},
               bios::LeafBioState{FT},
               k_ant::FT,
               k_brown::FT,
               k_cab::FT,
               k_car_v::FT,
               k_car_z::FT,
               k_cbc::FT,
               k_H₂O::FT,
               k_lma::FT,
               k_pro::FT,
               lwc::FT,
               x::FT,
               N::Int) where {FT}

Return the fraction of chlorophyll a and b absorption in the sublayer, given
- `bios` leaf biophysical state variables
- `k_ant` anthocyanin absorption coefficient
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
- `N` number of sublayers of each layer

"""
function sublayer_τ(
            biot::LeafBioTrait{FT},
            bios::LeafBioState{FT},
            k_ant::FT,
            k_brown::FT,
            k_cab::FT,
            k_car_v::FT,
            k_car_z::FT,
            k_cbc::FT,
            k_H₂O::FT,
            k_lma::FT,
            k_pro::FT,
            lwc::FT,
            x::FT,
            N::Int) where {FT}
    Σkx = k_ant * biot.ant +
          k_brown * biot.brown +
          k_cab * biot.cab +
          k_car_v * biot.car * (1 - bios.f_zeax) +
          k_car_z * biot.car * bios.f_zeax +
          k_cbc * biot.cbc +
          k_H₂O * (lwc * M_H₂O(FT) / ρ_H₂O(FT) * 100) +
          k_lma * (biot.lma - biot.cbc - biot.pro) +
          k_pro * biot.pro;
    Σkx *= x;
    τ_layer = (1 - Σkx) * exp(-Σkx) + Σkx^2 * expint(Σkx + eps(FT));

    return τ_layer ^ (FT(1) / N)
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-18: add function to compute how much PPAR is absorbed by PSII (used for SIF and PPAR)
#
#######################################################################################################################################################################################################
"""

    psii_fraction(wl::FT; f_max::FT = FT(0.7)) where {FT}

Return the fraction of PSII PPAR absorption, given
- `wl` excitation wavelength

Note if you want to customize the contribution from PSII, you can overwrite this function externally.

"""
function psii_fraction(wl::FT; f_max::FT = FT(0.7)) where {FT}
    if wl < 670
        return f_max
    else
        return max(0, f_max - (wl - 670) / 80 * f_max);
    end;
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to update the sublayer absorption and transmittance within LeafBio
#     2023-Sep-15: save f_sife as well as f_cab and f_car
#     2023-Sep-16: fix a typo related to f_car (was using 1 - f_car)
#     2023-Sep-18: save τ_all_i after computing τ_sub_i
#     2023-Sep-19: save f_ppar at the same time
#
#######################################################################################################################################################################################################
"""

    leaf_sublayer_f_τ!(config::SPACConfiguration{FT}, bio::LeafBio{FT}, lwc::FT, N::Int) where {FT}
    leaf_sublayer_f_τ!(spectra::ReferenceSpectra{FT}, bio::LeafBio{FT}, lwc::FT, N::Int) where {FT}

Update the sublayer absorption and transmittance within `bio`, given
- `config` SPAC configuration
- `bio` LeafBio struct
- `lwc` leaf water content
- `N` number of sublayers of each layer
- `spectra` ReferenceSpectra struct

"""
function leaf_sublayer_f_τ!(config::SPACConfiguration{FT}, bio::LeafBio{FT}, lwc::FT, N::Int) where {FT}
    (; K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, Λ) = config.SPECTRA;

    x = 1 / bio.trait.meso_n;
    bio.auxil.f_cab  .= sublayer_f_cab.((bio.trait,), (bio.state,), K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, lwc);
    bio.auxil.f_car  .= sublayer_f_car.((bio.trait,), (bio.state,), K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, lwc);
    bio.auxil.f_ppar .= bio.auxil.f_cab .+ bio.auxil.f_car .* bio.state.ϕ_car_ppar;
    bio.auxil.f_psii .= psii_fraction.(Λ);
    bio.auxil.f_sife .= bio.auxil.f_cab .+ bio.auxil.f_car .* bio.state.ϕ_car;

    bio.auxil.τ_sub_1 .= sublayer_τ.((bio.trait,), (bio.state,), K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, lwc, x, N);
    bio.auxil.τ_sub_2 .= sublayer_τ.((bio.trait,), (bio.state,), K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, lwc, 1-x, N);
    bio.auxil.τ_all_1 .= bio.auxil.τ_sub_1 .^ N;
    bio.auxil.τ_all_2 .= bio.auxil.τ_sub_2 .^ N;

    return nothing
end;
