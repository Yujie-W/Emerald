#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-16: add function to update the SIF conversion matrix of the first layer
#     2023-Sep-16: clear vec_b and vec_f before the calculation (for the case of recalculating the SIF conversion matrix)
#
#######################################################################################################################################################################################################
"""

    layer_1_sif_vec!(τ_i_θ::FT, τ_i_12::FT, τ_i_21::FT, τ_sub::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, f_sife::FT, τ_sub_sif::SubArray, vec_b::SubArray, vec_f::SubArray, ϕ::SubArray, N::Int) where {FT}

Update the SIF conversion matrix of the first layer, given
- `τ_i_θ` transmittance of the incoming radiation at the air-water interface
- `τ_i_12` transmittance of the isotropic radiation at the air-water interface
- `τ_i_21` transmittance of the isotropic radiation at the water-air interface
- `τ_sub` transmittance within a sublayer
- `τ_θ` transmittance of the incoming radiation across the leaf layer
- `ρ_1` reflectance of the first layer (1)
- `ρ_2` reflectance of the second layer (n-1)
- `f_sife` SIF excitation scaling factor (f_cab + f_car * ϕ_car)
- `τ_sub_sif` transmittance within a sublayer for SIF
- `vec_b` SIF vector backwards (SubArray of a Matrix)
- `vec_f` SIF vector forwards (SubArray of a Matrix)
- `ϕ` SIF emission PDF (SubArray of a Vector)
- `N` number of sublayers of each layer

"""
function layer_1_sif_vec!(τ_i_θ::FT, τ_i_12::FT, τ_i_21::FT, τ_sub::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, f_sife::FT, τ_sub_sif::SubArray, vec_b::SubArray, vec_f::SubArray, ϕ::SubArray, N::Int) where {FT}
    # parameters required for the calculation that can be derived from the input parameters
    ρ_i_21 = 1 - τ_i_21;
    α_sub = 1 - τ_sub;
    τ_all = τ_sub ^ N;

    # clear the SIF vector
    vec_b .= 0;
    vec_f .= 0;

    # 1. note that we do redo the calculation of the reflectance and transmittance within the leaf layer so that we can better model the sif emission and escape
    #    here the radiation is directly from the environment, as the reflection from the first layer will not hit back, no scaling of the light source is required
    #    however, the rescaling from the reflectance within the layer is required (we will do it at the last step)
    #    for the emitted SIF, we need to account for the reabsorption within the leaf layer
    #        - the SIF up will be downscaled by _t_sub ^ (_i - 0.5) times
    #        - the SIF down will be downscaled by _t_sub ^ (N - _i + 0.5) times
    #    ϕ_sif to account for the contribution to SIF excitation
    #    now the radiation goes from up to down
    rad_i = τ_i_θ * f_sife;
    for i in 1:N
        vec_b .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 2. then the radiation goes from down to up after hitting the water-air interface
    rad_i *= ρ_i_21;
    for i in N:-1:1
        vec_b .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 3. we also need to acount for the radiation reflected by the lower n-1 layer
    #    as this radiation is reflected back and forth, here we use the total radiation that goes in the direction into the first layer
    #        τ_θ * ρ_2 / (1 - ρ_1 * ρ_2)
    #    after accounting for the reflection and transmission within the layer
    #        the total radiation that transmits through is τ_θ * ρ_2 / (1 - ρ_1 * ρ_2) * τ_1
    #        the total radiation that reflects back is τ_θ * ρ_2 / (1 - ρ_1 * ρ_2) * ρ_1
    #    in this case, the total radiation is conserved
    #    now the radiation goes from down to up
    rad_i = τ_θ * ρ_2 / (1 - ρ_1 * ρ_2) * τ_i_12 * f_sife;
    for i in N:-1:1
        vec_b .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 4. then the radiation goes from up to down
    rad_i *= ρ_i_21;
    for i in 1:N
        vec_b .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 5. rescale it by 1 / (1 - ρ_21 * t_all * ρ_21 * t_all) for the SIFE
    denom = 1 - τ_all * ρ_i_21 * τ_all * ρ_i_21;
    vec_b ./= denom;
    vec_f ./= denom;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-16: add function to update the SIF conversion matrix of the n-1 layer
#
#######################################################################################################################################################################################################
"""

    layer_2_sif_vec!(τ_sub::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, τ_2::FT, f_sife::FT, τ_sub_sif::SubArray, vec_b::SubArray, vec_f::SubArray, ϕ::SubArray, N::Int) where {FT}

Update the SIF conversion matrix of the n-1 layer, given
- `τ_sub` transmittance within a sublayer
- `τ_θ` transmittance of the incoming radiation across the leaf layer
- `ρ_1` reflectance of the first layer (1)
- `ρ_2` reflectance of the second layer (n-1)
- `τ_2` transmittance of the second layer (n-1)
- `f_sife` SIF excitation scaling factor (f_cab + f_car * ϕ_car)
- `τ_sub_sif` transmittance within a sublayer for SIF
- `vec_b` SIF vector backwards (SubArray of a Matrix)
- `vec_f` SIF vector forwards (SubArray of a Matrix)
- `ϕ` SIF emission PDF (SubArray of a Vector)
- `N` number of sublayers of each layer

"""
function layer_2_sif_vec!(τ_sub::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, τ_2::FT, f_sife::FT, τ_sub_sif::SubArray, vec_b::SubArray, vec_f::SubArray, ϕ::SubArray, N::Int) where {FT}
    # parameters required for the calculation
    α_sub = 1 - τ_sub;
    τ_all = τ_sub ^ N;

    # clear the SIF vector
    vec_b .= 0;
    vec_f .= 0;

    # 1. here we consider the n-1 layers as one single layer, and the SIF transmission within this effective layer is same as the computed τ_sub_2
    #    then we need to rescale the interface ρ and τ for the effective layer so that the computed layer level ρ and τ are same as computed
    ρ_i_12 = effective_ρ_12(ρ_2, τ_2, τ_all);
    ρ_i_21 = effective_ρ_21(ρ_2, τ_2, τ_all);
    τ_i_12 = 1 - ρ_i_12;

    # 2. besides the the direction transmission from the first layer, we also need to account for the radiation reflected by the upper layer (the denominator 1 - ρ_1 * ρ_2)
    #    after accounting for the reflection and transmission within the layer
    #        the total radiation that transmits through is τ_θ / (1 - ρ_1 * ρ_2) * τ_2
    #        the total radiation that reflects back is τ_θ / (1 - ρ_1 * ρ_2) * ρ_1
    #    combine with that of the first layer, the total radiation is conserved
    #        the radiation that goes into the 1st layer is τ_θ * ρ_2 / (1 - ρ_1 * ρ_2) * (1 - ρ_1)
    #        the radiation that goes into the n-1 layer is τ_θ / (1 - ρ_1 * ρ_2) * (1 - ρ_2)
    #        the sum is then τ_θ (conserved)
    #    now the radiation goes from up to down
    rad_i = τ_θ / (1 - ρ_1 * ρ_2) * τ_i_12 * f_sife;
    for i in 1:N
        vec_b .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 3. then the radiation goes from down to up
    rad_i *= ρ_i_21;
    for i in N:-1:1
        vec_b .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub .* ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 4. rescale it by 1 / (1 - ρ_21 * t_all * ρ_21 * t_all) for the SIFE
    denom = 1 - τ_all * ρ_i_21 * τ_all * ρ_i_21;
    vec_b ./= denom;
    vec_f ./= denom;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-16: add function to compute the SIF emission backward
#
#######################################################################################################################################################################################################
"""

    leaf_sif_b(sif_b_1::FT, sif_f_1::FT, sif_b_2::FT, ρ_1::FT, τ_1::FT, ρ_2::FT) where {FT}

Return the SIF emission backward, given
- `sif_b_1` SIF emission backward from the first layer
- `sif_f_1` SIF emission forward from the first layer
- `sif_b_2` SIF emission backward from the second layer
- `ρ_1` reflectance of the first layer (1)
- `τ_1` transmittance of the first layer (1)
- `ρ_2` reflectance of the second layer (n-1)

"""
function leaf_sif_b(sif_b_1::FT, sif_f_1::FT, sif_b_2::FT, ρ_1::FT, τ_1::FT, ρ_2::FT) where {FT}
    # the sif that goes backward is the sum of
    #     - sif_b_1
    #     - sif_f_1 reflected by the 2nd layer and then transmit through the 1st layer (needs to be rescaled)
    #     - sif_b_2 that transmit through the n-1 layer (needs to be rescaled)
    denom = 1 - ρ_1 * ρ_2;

    return sif_b_1 + sif_f_1 * ρ_2 * τ_1 / denom + sif_b_2 * τ_1 / denom
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-16: add function to compute the SIF emission forward
#
#######################################################################################################################################################################################################
"""

    leaf_sif_f(sif_f_1::FT, sif_b_2::FT, sif_f_2::FT, ρ_1::FT, ρ_2::FT, τ_2::FT) where {FT}

Return the SIF emission forward, given
- `sif_f_1` SIF emission forward from the first layer
- `sif_b_2` SIF emission backward from the second layer
- `sif_f_2` SIF emission forward from the second layer
- `ρ_1` reflectance of the first layer (1)
- `ρ_2` reflectance of the second layer (n-1)
- `τ_2` transmittance of the second layer (n-1)

"""
function leaf_sif_f(sif_f_1::FT, sif_b_2::FT, sif_f_2::FT, ρ_1::FT, ρ_2::FT, τ_2::FT) where {FT}

    # the sif that goes forward is the sum of
    #     - sif_f_1 that transmit through the 2nd layer (needs to be rescaled)
    #     - sif_b_2 reflected by the 1st layer and then transmit through the 2nd layer (needs to be rescaled)
    #     - sif_f_2
    denom = 1 - ρ_1 * ρ_2;

    return sif_f_1 * τ_2 / denom + sif_b_2 * ρ_1 * τ_2 / denom + sif_f_2
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-16: add function to update the SIF conversion matrix of the leaf
#
#######################################################################################################################################################################################################
"""

    leaf_sif_matrices!(config::SPACConfiguration{FT}, bio::HyperLeafBio{FT}, N::Int) where {FT}
    leaf_sif_matrices!(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperLeafBio{FT}, N::Int) where {FT}

Update the SIF conversion matrix of the leaf, given
- `config` SPAC configuration
- `bio` leaf biophysics
- `N` number of sublayers of each layer
- `lha` leaf hyperspectral absorption
- `wls` wavelength set

"""
function leaf_sif_matrices! end;

leaf_sif_matrices!(config::SPACConfiguration{FT}, bio::HyperLeafBio{FT}, N::Int) where {FT} = leaf_sif_matrices!(config.LHA, config.WLSET, bio, N);

leaf_sif_matrices!(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperLeafBio{FT}, N::Int) where {FT} = (
    (; Φ_PS) = lha;
    (; IΛ_SIF, IΛ_SIFE) = wls;

    # update the SIF emission vector per excitation wavelength
    ϕ = view(Φ_PS, IΛ_SIF);
    ρ_1_sif = view(bio.auxil.ρ_layer_1, IΛ_SIF);
    τ_1_sif = view(bio.auxil.τ_layer_1, IΛ_SIF);
    ρ_2_sif = view(bio.auxil.ρ_layer_2, IΛ_SIF);
    τ_2_sif = view(bio.auxil.τ_layer_2, IΛ_SIF);
    τ_sub_sif_1 = view(bio.auxil.τ_sub_1, IΛ_SIF);
    τ_sub_sif_2 = view(bio.auxil.τ_sub_2, IΛ_SIF);
    for i in eachindex(IΛ_SIFE)
        ii = IΛ_SIFE[i];
        vec_b_1 = view(bio.auxil.mat_b_1, :, i);
        vec_f_1 = view(bio.auxil.mat_f_1, :, i);
        vec_b_2 = view(bio.auxil.mat_b_2, :, i);
        vec_f_2 = view(bio.auxil.mat_f_2, :, i);
        vec_b   = view(bio.auxil.mat_b, :, i);
        vec_f   = view(bio.auxil.mat_f, :, i);
        τ_i_θ   = bio.auxil.τ_interface_θ[ii];      # the transmittance of the incoming radiation at the air-water interface
        τ_i_12  = bio.auxil.τ_interface_12[ii];     # the transmittance of the isotropic radiation at the air-water interface
        τ_i_21  = bio.auxil.τ_interface_21[ii];     # the transmittance of the isotropic radiation at the water-air interface
        τ_sub_1 = bio.auxil.τ_sub_1[ii];            # the transmittance within a sublayer of layer 1
        τ_sub_2 = bio.auxil.τ_sub_2[ii];            # the transmittance within a sublayer of layer 2 (n-1)
        τ_l_θ   = bio.auxil.τ_layer_θ[ii];          # the transmittance of the incoming radiation across the leaf layer 1
        ρ_l_1   = bio.auxil.ρ_layer_1[ii];          # the reflectance of isotropic radiation across layer 1
        ρ_l_2   = bio.auxil.ρ_layer_2[ii];          # the reflectance of isotropic radiation across layer 2 (n-1)
        τ_l_2   = bio.auxil.τ_layer_2[ii];          # the transmittance of isotropic radiation across layer 2 (n-1)
        f_sife  = bio.auxil.f_sife[ii];

        # update the SIF conversion matrix of the two layers
        layer_1_sif_vec!(τ_i_θ, τ_i_12, τ_i_21, τ_sub_1, τ_l_θ, ρ_l_1, ρ_l_2, f_sife, τ_sub_sif_1, vec_b_1, vec_f_1, ϕ, N);
        layer_2_sif_vec!(τ_sub_2, τ_l_θ, ρ_l_1, ρ_l_2, τ_l_2, f_sife, τ_sub_sif_2, vec_b_2, vec_f_2, ϕ, N);

        # compute the SIF emission vector backward and forward
        vec_b .= leaf_sif_b.(vec_b_1, vec_f_1, vec_b_2, ρ_1_sif, τ_1_sif, ρ_2_sif);
        vec_f .= leaf_sif_f.(vec_f_1, vec_b_2, vec_f_2, ρ_1_sif, ρ_2_sif, τ_2_sif);
    end;

    return nothing
);
