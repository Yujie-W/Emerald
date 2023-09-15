
# function to compute reflectance and transmittance of the single layer
#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute reflectance and transmittance of the single layer of a leaf
#
#######################################################################################################################################################################################################
"""

    layer_ρ_τ(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, x::FT, θ::FT; N::Int = 10) where {FT}
    layer_ρ_τ(lha::HyperspectralAbsorption{FT}, θ::FT) where {FT}

Return the reflectance and transmittance of the single layer of a leaf (if bio, lwc, x are not provided, return those for the case of no absorption), given
- `lha` leaf hyperspectral absorption coefficients
- `bio` leaf hyperspectral biophysics
- `lwc` leaf water content
- `x` proportion of the layer of the whole leaf
- `θ` average incident angle of the incoming radiation
- `N` number of sublayers of each layer

"""
function layer_ρ_τ end;

layer_ρ_τ(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, x::FT, θ::FT; N::Int = 10) where {FT} = (
    (; NR) = lha;

    # radiation beam from the surface (100%)
    _τ_surface = interface_isotropic_τ.(FT(1), NR, θ);
    _ρ_surface = 1 .- _τ_surface;

    # scatter light from water to air (100%)
    _τ_21 = interface_isotropic_τ.(NR, FT(1), FT(90));
    _ρ_21 = 1 .- _τ_21;

    # absorption ratio within the layer
    # this one should be similar to sublayer_transmittance(lha, bio, x, 1)
    _t_all = sublayer_τ(lha, bio, lwc, x, N)[1] .^ N;

    # compute the integrated reflectance
    _denom = 1 .- _ρ_21 .* _t_all .* _ρ_21 .* _t_all;
    _ρ_int = _ρ_surface .+ _τ_surface .* _t_all .* _ρ_21 .* _t_all .* _τ_21 ./ _denom;
    _τ_int = _τ_surface .* _t_all .* _τ_21 ./ _denom;

    # You can double check the conservation like this
    # _α_int = _τ_surface .* (1 .- _t_all) ./ (1 .- _t_all .* _ρ_21);

    return _ρ_int, _τ_int
);

layer_ρ_τ(lha::HyperspectralAbsorption{FT}, θ::FT) where {FT} = (
    (; NR) = lha;

    # radiation beam from the surface (100%)
    _τ_surface = interface_isotropic_τ.(FT(1), NR, θ);
    _ρ_surface = 1 .- _τ_surface;

    # scatter light from water to air (100%)
    _τ_21 = interface_isotropic_τ.(NR, FT(1), FT(90));
    _ρ_21 = 1 .- _τ_21;

    # absorption ratio within the layer
    # this one should be similar to sublayer_transmittance(lha, bio, x, 1)
    _t_all = ones(length(NR));

    # compute the integrated reflectance
    _denom = 1 .- _ρ_21 .* _t_all .* _ρ_21 .* _t_all;
    _ρ_int = _ρ_surface .+ _τ_surface .* _t_all .* _ρ_21 .* _t_all .* _τ_21 ./ _denom;
    _τ_int = _τ_surface .* _t_all .* _τ_21 ./ _denom;

    # You can double check the conservation like this
    # _α_int = _τ_surface .* (1 .- _t_all) ./ (1 .- _t_all .* _ρ_21);

    return _ρ_int, _τ_int
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute reflectance and transmittance of the single layer of a leaf with non isotropic radiation
#
#######################################################################################################################################################################################################
"""

    layer_ρ_τ_direct(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT; N::Int = 10) where {FT}

Return the reflectance and transmittance of the single layer of a leaf with non isotropic radiation, given
- `lha` leaf hyperspectral absorption coefficients
- `bio` leaf hyperspectral biophysics
- `lwc` leaf water content
- `θ` average incident angle of the incoming radiation
- `N` number of sublayers of each layer

"""
function layer_ρ_τ_direct(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT; N::Int = 10) where {FT}
    (; NR) = lha;
    (; MESOPHYLL_N) = bio;

    # radiation beam from the surface (100%)
    _τ_surface = interface_isotropic_τ.(FT(1), NR, θ);
    _ρ_surface = 1 .- _τ_surface;

    # scatter light from water to air (100%)
    _τ_21 = interface_isotropic_τ.(NR, FT(1), FT(90));
    _ρ_21 = 1 .- _τ_21;

    # absorption ratio within the layer
    # this one should be similar to sublayer_transmittance(lha, bio, x, 1)
    _t_all = sublayer_τ(lha, bio, lwc, 1/MESOPHYLL_N, N)[1] .^ N;

    # compute the integrated reflectance
    _denom = 1 .- _ρ_21 .* _t_all .* _ρ_21 .* _t_all;
    _ρ_int = _ρ_surface .+ _τ_surface .* _t_all .* _ρ_21 .* _t_all .* _τ_21 ./ _denom;
    _τ_int = _τ_surface .* _t_all .* _τ_21 ./ _denom;

    # You can double check the conservation like this
    # _α_int = _τ_surface .* (1 .- _t_all) ./ (1 .- _t_all .* _ρ_21);

    return _ρ_int, _τ_int
end;

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

    leaf_layer_spectra!(bio::HyperLeafBio{FT}, N::Int) where {FT}

Update the reflectance and transmittance of the leaf layers, given
- `bio` leaf hyperspectral biophysics
- `N` number of sublayers of each layer

"""
function leaf_layer_spectra!(bio::HyperLeafBio{FT}, N::Int) where {FT}
    bio.auxil.ρ_layer_θ .= layer_1_ρ.(bio.auxil.τ_interface_θ , bio.auxil.ρ_21, bio.auxil.τ_sub_1, N);
    bio.auxil.τ_layer_θ .= layer_1_τ.(bio.auxil.τ_interface_θ , bio.auxil.ρ_21, bio.auxil.τ_sub_1, N);
    bio.auxil.ρ_layer_1 .= layer_1_ρ.(bio.auxil.τ_interface_12, bio.auxil.ρ_21, bio.auxil.τ_sub_1, N);
    bio.auxil.τ_layer_1 .= layer_1_τ.(bio.auxil.τ_interface_12, bio.auxil.ρ_21, bio.auxil.τ_sub_1, N);

    x = 1 - 1 / bio.state.meso_n;
    bio.auxil.ρ_layer_2 .= layer_2_ρ.(bio.auxil.ρ_layer_1, bio.auxil.τ_layer_1, x, N);
    bio.auxil.τ_layer_2 .= layer_2_τ.(bio.auxil.ρ_layer_1, bio.auxil.τ_layer_1, x, N);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute reflectance and transmittance of the single layer of a leaf with diffuse radiation
#
#######################################################################################################################################################################################################
"""

    layer_ρ_τ_diffuse(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT; N::Int = 10) where {FT}
    layer_ρ_τ_diffuse(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}) where {FT}

Return the reflectance and transmittance of the single layer of a leaf with diffuse radiation (if lwc and N not provided, return the values of the case with zero absorption), given
- `lha` leaf hyperspectral absorption coefficients
- `bio` leaf hyperspectral biophysics
- `lwc` leaf water content
- `N` number of sublayers of each layer

"""
function layer_ρ_τ_diffuse end;

layer_ρ_τ_diffuse(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT; N::Int = 10) where {FT} = (
    (; NR) = lha;
    (; MESOPHYLL_N) = bio;

    # radiation beam from the surface (100%)
    _τ_surface = interface_isotropic_τ.(FT(1), NR, FT(90));
    _ρ_surface = 1 .- _τ_surface;

    # scatter light from water to air (100%)
    _τ_21 = interface_isotropic_τ.(NR, FT(1), FT(90));
    _ρ_21 = 1 .- _τ_21;

    # absorption ratio within the layer
    # this one should be similar to sublayer_transmittance(lha, bio, x, 1)
    _t_all = sublayer_τ(lha, bio, lwc, 1/MESOPHYLL_N, N)[1] .^ N;

    # compute the integrated reflectance for the first layer
    _ρ_1 = similar(_t_all);
    _τ_1 = similar(_t_all);
    _ρ_2 = similar(_t_all);
    _τ_2 = similar(_t_all);
    _denom = similar(_t_all);
    @. _denom = 1 - _ρ_21 * _t_all * _ρ_21 * _t_all;
    @. _ρ_1 = _ρ_surface + _τ_surface * _t_all * _ρ_21 * _t_all * _τ_21 / _denom;
    @. _τ_1 = _τ_surface * _t_all * _τ_21 / _denom;

    # compute the integrated reflectance and transmittance for the n-1 layer
    d = similar(_t_all);
    a = similar(_t_all);
    b = similar(_t_all);
    @. d = sqrt((_τ_1 ^ 2 - _ρ_1 ^ 2 - 1) ^ 2 - 4 * _ρ_1 ^ 2);
    @. a = (1 + _ρ_1 ^ 2 - _τ_1 ^ 2 + d) / (2 * _ρ_1);
    @. b = (1 - _ρ_1 ^ 2 + _τ_1 ^ 2 + d) / (2 * _τ_1);
    m = MESOPHYLL_N - 1;
    @. _ρ_2 = (b ^ m - 1 / (b ^ m)) / (a * b ^ m - 1 / (a * b ^ m));
    @. _τ_2 = (a - 1 / a) / (a * b ^ m - 1 / (a * b ^ m));

    return _ρ_1, _τ_1, _ρ_2, _τ_2
);

layer_ρ_τ_diffuse(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}) where {FT} = (
    (; NR) = lha;
    (; MESOPHYLL_N) = bio;

    # radiation beam from the surface (100%)
    _τ_surface = interface_isotropic_τ.(FT(1), NR, FT(90));
    _ρ_surface = 1 .- _τ_surface;

    # scatter light from water to air (100%)
    _τ_21 = interface_isotropic_τ.(NR, FT(1), FT(90));
    _ρ_21 = 1 .- _τ_21;

    # absorption ratio within the layer
    # this one should be similar to sublayer_transmittance(lha, bio, x, 1)
    _t_all = ones(length(NR));

    # compute the integrated reflectance for the first layer
    _ρ_1 = similar(_t_all);
    _τ_1 = similar(_t_all);
    _ρ_2 = similar(_t_all);
    _τ_2 = similar(_t_all);
    _denom = similar(_t_all);
    @. _denom = 1 - _ρ_21 * _t_all * _ρ_21 * _t_all;
    @. _ρ_1 = _ρ_surface + _τ_surface * _t_all * _ρ_21 * _t_all * _τ_21 / _denom;
    @. _τ_1 = _τ_surface * _t_all * _τ_21 / _denom;

    # compute the integrated reflectance and transmittance for the n-1 layer (d = 0, a = 1, b = 1 here after solved analytically)
    #     solve this:
    #         tau(n) = tau(1) * tau(n-1) / (1 - rho(1) * rho(n-1))                  =>
    #         tau(n) = tau(1) * tau(n-1) / (tau(1) + tau(n-1) - tau(1) * tau(n-1))  =>
    # according to wolframalpha.com
    #         tau(n) = tau(1) / (tau(1) + n * tau(1))
    #
    m = MESOPHYLL_N - 1;
    @. _τ_2 = _τ_1 / (_τ_1 + m - m * _τ_1);
    @. _ρ_2 = 1 - _τ_2;

    return _ρ_1, _τ_1, _ρ_2, _τ_2
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute the SIF excitation ratio of a single layer of a leaf
# To do
#     TODO: this function is not used, consider to remove it later
#
#######################################################################################################################################################################################################
"""

    layer_raw_sife(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, x::FT, θ::FT; N::Int = 10) where {FT}

Return the SIF excitation ratio of a single layer of a leaf, given
- `lha` leaf hyperspectral absorption coefficients
- `wls` hyperspectral wavelength set
- `bio` leaf hyperspectral biophysics
- `lwc` leaf water content
- `x` proportion of the layer of the whole leaf
- `θ` average incident angle of the incoming radiation
- `N` number of sublayers of each layer

"""
function layer_raw_sife(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, x::FT, θ::FT; N::Int = 10) where {FT}
    (; NR) = lha;
    (; IΛ_SIFE) = wls;

    # radiation beam from the surface (100%)
    _τ_surface = interface_isotropic_τ.(FT(1), NR[IΛ_SIFE], θ);

    # scatter light from water to air (100%)
    _τ_21 = interface_isotropic_τ.(NR[IΛ_SIFE], FT(1), FT(90));
    _ρ_21 = 1 .- _τ_21;

    # compute the transmittance within the sublayer
    _rad_i = deepcopy(_τ_surface);
    _t_sub = sublayer_τ(lha, bio, lwc, x, N)[1][IΛ_SIFE];
    _a_sub = 1 .- _t_sub;
    _t_all = _t_sub .^ N;

    # itereate through the sublayers and add up the SIF excitation ratio (TODO: consinder reabsorption later)
    _sif_excitation = zeros(FT, length(IΛ_SIFE));
    for _ in 1:N
        _sif_excitation += _rad_i .* _a_sub;
        _rad_i .*= _t_sub;
    end;

    # rescale the matrices to account for the 1-2 and 2-1 reflection cycle
    _denom = 1 .- _ρ_21 .* _t_all;
    _sif_excitation ./= _denom;

    return _sif_excitation
end;
