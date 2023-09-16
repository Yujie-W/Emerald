#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the rescaled reflectance of the first layer of a leaf
#
#######################################################################################################################################################################################################
"""

    layer_1_ρ(τ_in::FT, ρ_21::FT, τ_sub::FT, N::Int) where {FT}
    layer_1_ρ(τ_in::FT, ρ_21::FT, τ_all::FT) where {FT}

Return the rescaled reflectance of the first layer of a leaf, given
- `τ_in` transmittance of the incoming radiation
- `ρ_21` reflectance at the water(2)-air(1) interface
- `τ_sub` transmittance within a sublayer
- `N` number of sublayers of each layer
- `τ_all` transmittance within a layer

"""
function layer_1_ρ end;

layer_1_ρ(τ_in::FT, ρ_21::FT, τ_sub::FT, N::Int) where {FT} = layer_1_ρ(τ_in, ρ_21, τ_sub ^ N);

layer_1_ρ(τ_in::FT, ρ_21::FT, τ_all::FT) where {FT} = (
    denom = 1 - ρ_21 * τ_all * ρ_21 * τ_all;

    return 1 - τ_in + τ_in * τ_all * ρ_21 * τ_all * (1 - ρ_21) / denom
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the rescaled transmittance of the first layer of a leaf
#
#######################################################################################################################################################################################################
"""

    layer_1_τ(τ_in::FT, ρ_21::FT, τ_sub::FT, N::Int) where {FT}
    layer_1_τ(τ_in::FT, ρ_21::FT, τ_all::FT) where {FT}

Return the rescaled transmittance of the first layer of a leaf, given
- `τ_in` transmittance of the incoming radiation
- `ρ_21` reflectance at the water(2)-air(1) interface
- `τ_sub` transmittance within a sublayer
- `N` number of sublayers of each layer
- `τ_all` transmittance within a layer

"""
function layer_1_τ end;

layer_1_τ(τ_in::FT, ρ_21::FT, τ_sub::FT, N::Int) where {FT} = layer_1_τ(τ_in, ρ_21, τ_sub ^ N);

layer_1_τ(τ_in::FT, ρ_21::FT, τ_all::FT) where {FT} = (
    denom = 1 - ρ_21 * τ_all * ρ_21 * τ_all;

    return τ_in * τ_all * (1 - ρ_21) / denom
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the rescaled reflectance of the n-1 layer of a leaf
#
# Note, when there is no absorption, d = 0, a = 1, and b = 1 here.
# Then the problem to compute the integrated reflectance and transmittance for the n-1 layer is to solve
#     tau(n) = tau(1) * tau(n-1) / (1 - rho(1) * rho(n-1))                  =>
#     tau(n) = tau(1) * tau(n-1) / (tau(1) + tau(n-1) - tau(1) * tau(n-1))  =>
# according to wolframalpha.com
#      tau(n) = tau(1) / (tau(1) + n * tau(1))
#
#######################################################################################################################################################################################################
"""

    layer_2_ρ(ρ₁::FT, τ₁::FT, m::FT) where {FT}

Return the rescaled reflectance of the n-1 layer of a leaf, given
- `ρ₁` reflectance of the first layer
- `τ₁` transmittance of the first layer
- `m` n - 1

"""
function layer_2_ρ(ρ₁::FT, τ₁::FT, m::FT) where {FT}
    ρ₁² = ρ₁ ^ 2;
    τ₁² = τ₁ ^ 2;
    d = sqrt((τ₁² - ρ₁² - 1) ^ 2 - 4ρ₁²);
    a = (1 + ρ₁² - τ₁² + d) / 2ρ₁;
    b = (1 - ρ₁² + τ₁² + d) / 2τ₁;

    return (b ^ m - 1 / (b ^ m)) / (a * b ^ m - 1 / (a * b ^ m))
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the rescaled transmittance of the n-1 layer of a leaf
#
#######################################################################################################################################################################################################
"""

    layer_2_τ(ρ₁::FT, τ₁::FT, m::FT) where {FT}

Return the rescaled reflectance of the n-1 layer of a leaf, given
- `ρ₁` reflectance of the first layer
- `τ₁` transmittance of the first layer
- `m` n - 1

"""
function layer_2_τ(ρ₁::FT, τ₁::FT, m::FT) where {FT}
    ρ₁² = ρ₁ ^ 2;
    τ₁² = τ₁ ^ 2;
    d = sqrt((τ₁² - ρ₁² - 1) ^ 2 - 4ρ₁²);
    a = (1 + ρ₁² - τ₁² + d) / 2ρ₁;
    b = (1 - ρ₁² + τ₁² + d) / 2τ₁;

    return (a - 1 / a) / (a * b ^ m - 1 / (a * b ^ m))
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to update the reflectance and transmittance of the leaf layers
#
#######################################################################################################################################################################################################
"""

    leaf_layer_ρ_τ!(bio::HyperLeafBio{FT}, N::Int) where {FT}

Update the reflectance and transmittance of the leaf layers, given
- `bio` leaf hyperspectral biophysics
- `N` number of sublayers of each layer

"""
function leaf_layer_ρ_τ!(bio::HyperLeafBio{FT}, N::Int) where {FT}
    bio.auxil.ρ_layer_θ .= layer_1_ρ.(bio.auxil.τ_interface_θ , bio.auxil.ρ_interface_21, bio.auxil.τ_sub_1, N);
    bio.auxil.τ_layer_θ .= layer_1_τ.(bio.auxil.τ_interface_θ , bio.auxil.ρ_interface_21, bio.auxil.τ_sub_1, N);
    bio.auxil.ρ_layer_1 .= layer_1_ρ.(bio.auxil.τ_interface_12, bio.auxil.ρ_interface_21, bio.auxil.τ_sub_1, N);
    bio.auxil.τ_layer_1 .= layer_1_τ.(bio.auxil.τ_interface_12, bio.auxil.ρ_interface_21, bio.auxil.τ_sub_1, N);

    m = bio.state.meso_n - 1;
    bio.auxil.ρ_layer_2 .= layer_2_ρ.(bio.auxil.ρ_layer_1, bio.auxil.τ_layer_1, m);
    bio.auxil.τ_layer_2 .= layer_2_τ.(bio.auxil.ρ_layer_1, bio.auxil.τ_layer_1, m);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-16: add function to update the SIF conversion matrix of the first layer
#
#######################################################################################################################################################################################################
"""

    layer_1_sif_vec!(
                τ_i_θ::FT,
                τ_i_12::FT,
                τ_i_21::FT,
                τ_sub::FT,
                τ_θ::FT,
                ρ_1::FT,
                ρ_2::FT,
                f_sife::FT,
                τ_sub_sif::SubArray,
                mat_b::SubArray,
                mat_f::SubArray,
                Φ_PS::SubArray,
                N::Int) where {FT}

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
- `Φ_PS` SIF emission PDF (SubArray of a Vector)
- `N` number of sublayers of each layer

"""
function layer_1_sif_vec!(
            τ_i_θ::FT,
            τ_i_12::FT,
            τ_i_21::FT,
            τ_sub::FT,
            τ_θ::FT,
            ρ_1::FT,
            ρ_2::FT,
            f_sife::FT,
            τ_sub_sif::SubArray,
            vec_b::SubArray,
            vec_f::SubArray,
            Φ_PS::SubArray,
            N::Int) where {FT}
    # parameters required for the calculation that can be derived from the input parameters
    ρ_i_21 = 1 - τ_i_21;
    α_sub = 1 - τ_sub;
    τ_all = τ_sub ^ N;

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
        vec_b .+= rad_i / 2 * α_sub * Φ_PS .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub * Φ_PS .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 2. then the radiation goes from down to up after hitting the water-air interface
    rad_i *= ρ_i_21;
    for i in N:-1:1
        vec_b .+= rad_i / 2 * α_sub * Φ_PS .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub * Φ_PS .* τ_sub_sif .^ (N - i + FT(0.5));
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
        vec_b .+= rad_i / 2 * α_sub * Φ_PS .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub * Φ_PS .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 4. then the radiation goes from up to down
    rad_i *= ρ_i_21;
    for i in 1:N
        vec_b .+= rad_i / 2 * α_sub * Φ_PS .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub * Φ_PS .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 5. rescale it by 1 / (1 - ρ_21 * t_all * ρ_21 * t_all) for the SIFE
    denom = 1 - τ_all * ρ_i_21 * τ_all * ρ_i_21;
    vec_b ./= denom;
    vec_f ./= denom;

    return nothing
end;
