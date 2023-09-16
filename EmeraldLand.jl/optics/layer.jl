#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the rescaled reflectance of the first layer of a leaf
#
#######################################################################################################################################################################################################
"""

    layer_1_ρ(τ_in::FT, ρ_21::FT, τ_sub::FT, N::Int) where {FT}

Return the rescaled reflectance of the first layer of a leaf, given
- `τ_in` transmittance of the incoming radiation
- `ρ_21` reflectance at the water(2)-air(1) interface
- `τ_sub` transmittance within a sublayer
- `N` number of sublayers of each layer

"""
function layer_1_ρ(τ_in::FT, ρ_21::FT, τ_sub::FT, N::Int) where {FT}
    τ_all = τ_sub ^ N;
    denom = 1 - ρ_21 * τ_all * ρ_21 * τ_all;

    return 1 - τ_in + τ_in * τ_all * ρ_21 * τ_all * (1 - ρ_21) / denom
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the rescaled transmittance of the first layer of a leaf
#
#######################################################################################################################################################################################################
"""

    layer_1_τ(τ_in::FT, ρ_21::FT, τ_sub::FT, N::Int) where {FT}

Return the rescaled transmittance of the first layer of a leaf, given
- `τ_in` transmittance of the incoming radiation
- `ρ_21` reflectance at the water(2)-air(1) interface
- `τ_sub` transmittance within a sublayer
- `N` number of sublayers of each layer

"""
function layer_1_τ(τ_in::FT, ρ_21::FT, τ_sub::FT, N::Int) where {FT}
    τ_all = τ_sub ^ N;
    denom = 1 - ρ_21 * τ_all * ρ_21 * τ_all;

    return τ_in * τ_all * (1 - ρ_21) / denom
end;


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
- `m` proportion of the n-1 layers of the whole leaf

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
- `m` proportion of the n-1 layers of the whole leaf

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
    bio.auxil.ρ_layer_θ .= layer_1_ρ.(bio.auxil.τ_interface_θ , bio.auxil.ρ_21, bio.auxil.τ_sub_1, N);
    bio.auxil.τ_layer_θ .= layer_1_τ.(bio.auxil.τ_interface_θ , bio.auxil.ρ_21, bio.auxil.τ_sub_1, N);
    bio.auxil.ρ_layer_1 .= layer_1_ρ.(bio.auxil.τ_interface_12, bio.auxil.ρ_21, bio.auxil.τ_sub_1, N);
    bio.auxil.τ_layer_1 .= layer_1_τ.(bio.auxil.τ_interface_12, bio.auxil.ρ_21, bio.auxil.τ_sub_1, N);

    x = 1 - 1 / bio.state.meso_n;
    bio.auxil.ρ_layer_2 .= layer_2_ρ.(bio.auxil.ρ_layer_1, bio.auxil.τ_layer_1, x, N);
    bio.auxil.τ_layer_2 .= layer_2_τ.(bio.auxil.ρ_layer_1, bio.auxil.τ_layer_1, x, N);

    return nothing
end;
