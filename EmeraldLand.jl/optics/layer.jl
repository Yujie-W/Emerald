
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
- `N` number of sublayers of the each sublayer (default: 10)

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
- `N` number of sublayers of the each sublayer (default: 10)

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
- `N` number of sublayers of the each sublayer (default: 10)

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
- `N` number of sublayers of the each layer (default: 10)

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
