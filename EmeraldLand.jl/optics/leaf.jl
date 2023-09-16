#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the reflectance of the leaf
#
#######################################################################################################################################################################################################
"""

    leaf_ρ(ρ₀::FT, τ₀::FT, ρ₁::FT, τ₁::FT, ρ₂::FT) where {FT}

Return the leaf level reflectance, given
- `ρ₀` reflectance of the first layer with incident radiation
- `τ₀` transmittance of the first layer with incident radiation
- `ρ₁` reflectance of the first layer with isotropic radiation
- `τ₁` transmittance of the first layer with isotropic radiation
- `ρ₂` reflectance of the n-1 layer with isotropic radiation

"""
function leaf_ρ(ρ₀::FT, τ₀::FT, ρ₁::FT, τ₁::FT, ρ₂::FT) where {FT}
    denom = 1 - ρ₁ * ρ₂;

    return ρ₀ + τ₀ * ρ₂ * τ₁ / denom
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the transmittance of the leaf
#
#######################################################################################################################################################################################################
"""

    leaf_τ(τ₀::FT, ρ₁::FT, ρ₂::FT, τ₂::FT) where {FT}

Return the leaf level reflectance, given
- `τ₀` transmittance of the first layer with incident radiation
- `ρ₁` reflectance of the first layer with isotropic radiation
- `ρ₂` reflectance of the n-1 layer with isotropic radiation
- `τ₂` transmittance of the n-1 layer with isotropic radiation

"""
function leaf_τ(τ₀::FT, ρ₁::FT, ρ₂::FT, τ₂::FT) where {FT}
    denom = 1 - ρ₁ * ρ₂;

    return τ₀ * τ₂ / denom
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to update the leaf reflectance and transmittance within HyperLeafBio
#
#######################################################################################################################################################################################################
"""

    leaf_ρ_τ!(bio::HyperLeafBio{FT}) where {FT}

Update the leaf reflectance and transmittance within `bio`, given
- `bio` HyperLeafBio struct

"""
function leaf_ρ_τ!(bio::HyperLeafBio{FT}) where {FT}
    bio.auxil.ρ_leaf .= leaf_ρ.(bio.auxil.ρ_layer_θ, bio.auxil.τ_layer_θ, bio.auxil.ρ_layer_1, bio.auxil.τ_layer_1, bio.auxil.ρ_layer_2);
    bio.auxil.τ_leaf .= leaf_τ.(bio.auxil.τ_layer_θ, bio.auxil.ρ_layer_1, bio.auxil.ρ_layer_2, bio.auxil.τ_layer_2);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to run all the step within one function all
#
#######################################################################################################################################################################################################
"""

    leaf_ρ_τ!(config::SPACConfiguration{FT}, bio::HyperLeafBio{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10) where {FT}

Update the interface, sublayer, layer, and leaf level reflectance and transmittance within `bio`, given
- `config` SPAC configuration
- `bio` HyperLeafBio struct

"""
function leaf_spectra! end;

leaf_spectra!(config::SPACConfiguration{FT}, bio::HyperLeafBio{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10) where {FT} = (
    leaf_interface_ρ_τ!(config, bio, θ);
    leaf_sublayer_f_τ!(config, bio, lwc, N);
    leaf_layer_ρ_τ!(bio, N);
    leaf_ρ_τ!(bio);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-16: add function to compute the leaf level SIF vectors
#
#######################################################################################################################################################################################################
"""

    leaf_sif_matrices!(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10, ϕ_car::FT = FT(0)) where {FT}

Return the leaf level SIF matrices (after leaf reabsorption), given
- `lha` leaf hyperspectral absorption coefficients
- `wls` hyperspectral wavelength set
- `bio` leaf hyperspectral biophysics
- `lwc` leaf water content
- `θ` incident angle of the incoming radiation
- `N` number of sublayers of each layer
- `ϕ_car` carotenoid contribution to chlorophyll fluorescence (default: 0)

"""
function leaf_sif_matrices! end;

leaf_sif_vector!(N::Int) where {FT} = (
    # parameters required for the calculation
    nr = 1.3;
    τ_i_θ = 0;
    τ_i_12 = 0;
    τ_i_21 = 0;
    ρ_i_21 = 1 - τ_i_21;
    τ_θ = 0;
    ρ_1 = 0;
    τ_1 = 0;
    ρ_2 = 0;
    τ_2 = 0;
    τ_sub_1 = 0;
    α_sub_1 = 1 - τ_sub_1;
    τ_all_1 = τ_sub_1 ^ N;
    τ_sub_2 = 0;
    α_sub_2 = 1 - τ_sub_2;
    τ_all_2 = τ_sub_2 ^ N;
    f_cab = 0;
    f_car = 0;
    ϕ_car = 0;
    ϕ_sif = f_cab + f_car * ϕ_car;

    # vector with the length of SIF wavelengths
    τ_sub_sif_1 = zeros(FT, length(IΛ_SIF));
    τ_sub_sif_2 = zeros(FT, length(IΛ_SIF));
    vec_b_1 = zeros(FT, length(IΛ_SIF));
    vec_b_2 = zeros(FT, length(IΛ_SIF));
    vec_f_1 = zeros(FT, length(IΛ_SIF));
    vec_f_2 = zeros(FT, length(IΛ_SIF));
    vec_b = zeros(FT, length(IΛ_SIF));
    vec_f = zeros(FT, length(IΛ_SIF));
    Φ_PS = zeros(FT, length(IΛ_SIF));
    ρ_1_sif = zeros(FT, length(IΛ_SIF));
    ρ_2_sif = zeros(FT, length(IΛ_SIF));
    τ_1_sif = zeros(FT, length(IΛ_SIF));
    τ_2_sif = zeros(FT, length(IΛ_SIF));
    denom_sif = zeros(FT, length(IΛ_SIF));

    #
    # PART 1: compute the matrices for the upper layer
    #

    # 1. now the radiation goes from up to down
    #    the SIF up will be downscaled by _t_sub ^ (_i - 0.5) times
    #    the SIF down will be downscaled by _t_sub ^ (N - _i + 0.5) times
    rad_i = τ_i_θ * ϕ_sif;
    for i in 1:N
        vec_b_1 .+= rad_i / 2 * α_sub_1 * Φ_PS .* τ_sub_sif_1 .^ (i - FT(0.5));
        vec_f_1 .+= rad_i / 2 * α_sub_1 * Φ_PS .* τ_sub_sif_1 .^ (N - i + FT(0.5));
        rad_i *= τ_sub_1;
    end;

    # 2. then the radiation goes from down to up
    rad_i *= ρ_i_21;
    for i in N:-1:1
        vec_b_1 .+= rad_i / 2 * α_sub_1 * Φ_PS .* τ_sub_sif_1 .^ (i - FT(0.5));
        vec_f_1 .+= rad_i / 2 * α_sub_1 * Φ_PS .* τ_sub_sif_1 .^ (N - i + FT(0.5));
        rad_i *= τ_sub_1;
    end;

    # 3. we also need to acount for the radiation reflected by the lower n-1 layer
    #    now the radiation goes from down to up
    rad_i = τ_θ * ρ_2 / (1 - ρ_1 * ρ_2) * τ_i_12 * ϕ_sif;
    for i in N:-1:1
        vec_b_1 .+= rad_i / 2 * α_sub_1 * Φ_PS .* τ_sub_sif_1 .^ (i - FT(0.5));
        vec_f_1 .+= rad_i / 2 * α_sub_1 * Φ_PS .* τ_sub_sif_1 .^ (N - i + FT(0.5));
        rad_i *= τ_sub_1;
    end;

    # 4. then the radiation goes from up to down
    rad_i *= ρ_i_21;
    for i in 1:N
        vec_b_1 .+= rad_i / 2 * α_sub_1 * Φ_PS .* τ_sub_sif_1 .^ (i - FT(0.5));
        vec_f_1 .+= rad_i / 2 * α_sub_1 * Φ_PS .* τ_sub_sif_1 .^ (N - i + FT(0.5));
        rad_i *= τ_sub_1;
    end;

    # 3. rescale it by 1 / (1 - ρ_21 * t_all * ρ_21 * t_all) for the SIFE
    denom = 1 - τ_all_1 * ρ_i_21 * τ_all_1 * ρ_i_21;
    vec_b_1 ./= denom;
    vec_f_1 ./= denom;


    #
    # Part 2: compute the matrices for the lower layer
    #

    # 1. here we consider the n-1 layers as one single layer, and the SIF transmission within this effective layer is same as the computed τ_sub_2
    #    then we need to rescale the interface ρ and τ for the effective layer so that the computed layer level ρ and τ are same as computed
    ρ_i_12_rescale = effective_ρ_12(ρ_2, τ_2, τ_all_2);
    ρ_i_21_rescale = effective_ρ_21(ρ_2, τ_2, τ_all_2);
    τ_i_12_rescale = 1 - ρ_i_12_rescale;

    # 2. compute the radiation from up to down
    #    we also need to account for the radiation reflected by the upper layer (the denominator)
    rad_i = τ_θ / (1 - ρ_1 * ρ_2) * τ_i_12_rescale * ϕ_sif;
    for i in 1:N
        vec_b_2 .+= rad_i / 2 * α_sub_2 * Φ_PS .* τ_sub_sif_2 .^ (i - FT(0.5));
        vec_f_2 .+= rad_i / 2 * α_sub_2 * Φ_PS .* τ_sub_sif_2 .^ (N - i + FT(0.5));
        rad_i *= τ_sub_2;
    end;

    # 3. then the radiation goes from down to up
    rad_i *= ρ_i_21_rescale;
    for i in N:-1:1
        vec_b_2 .+= rad_i / 2 * α_sub_2 * Φ_PS .* τ_sub_sif_2 .^ (i - FT(0.5));
        vec_f_2 .+= rad_i / 2 * α_sub_2 * Φ_PS .* τ_sub_sif_2 .^ (N - i + FT(0.5));
        rad_i *= τ_sub_2;
    end;

    # 4. rescale it by 1 / (1 - ρ_21 * t_all * ρ_21 * t_all) for the SIFE
    denom = 1 - τ_all_2 * ρ_i_21_rescale * τ_all_2 * ρ_i_21_rescale;
    vec_b_2 ./= denom;
    vec_f_2 ./= denom;

    #
    # Part 3: compute the matrices for the entire leaf
    #

    # in part 2, we already computed the SIF vectors that escape from the up and lower interfaces of the two layers (1 and n-1)
    # here we need to consider the reflectance of SIF between the two layers to compute the matrices for the entire leaf (forward and backward)
    # the sif that goes backward is the sum of
    #     - vec_b_1
    #     - vec_f_1 reflected by the 2nd layer and then transmit through the 1st layer (needs to be rescaled)
    #     - vec_b_2 that transmit through the n-1 layer (needs to be rescaled)
    # the sif that goes forward is the sum of
    #     - vec_f_1 that transmit through the 2nd layer (needs to be rescaled)
    #     - vec_b_2 reflected by the 1st layer and then transmit through the 2nd layer (needs to be rescaled)
    #     - vec_f_2
    denom_sif .= 1 .- ρ_1_sif .* ρ_2_sif;
    vec_b .= vec_b_1 .+ vec_f_1 .* ρ_2_sif .* τ_1_sif ./ denom .+ vec_b_2 .* τ_1_sif ./ denom;
    vec_f .= vec_f_1 .* τ_2_sif ./ denom .+ vec_b_2 .* ρ_1_sif .* τ_2_sif ./ denom .+ vec_f_2;

    return nothing
);
