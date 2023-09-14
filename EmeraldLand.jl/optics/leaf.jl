#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute the leaf spectra for reflectance and transmittance
#
#######################################################################################################################################################################################################
"""

    leaf_spectra(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10) where {FT}

Return the leaf spectra for reflectance and transmittance, given
- `lha` leaf hyperspectral absorption coefficients
- `bio` leaf hyperspectral biophysics
- `lwc` leaf water content
- `θ` incident angle of the incoming radiation
- `N` number of sublayers of the each sublayer (default: 10)

"""
function leaf_spectra(lha::HyperspectralAbsorption{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10) where {FT}
    (; MESOPHYLL_N) = bio;

    # compute the reflectance and transmittance of the adaxial and abaixal surfaces (with optimum angle)
    _ρ_s,_τ_s = layer_ρ_τ(lha, bio, lwc, 1/MESOPHYLL_N, θ; N= N);

    # define the vectors
    _denom = similar(_ρ_s);
    _ρ_leaf = similar(_ρ_s);
    _τ_leaf = similar(_ρ_s);

    # compute the reflectance and transmittance of the leaf layers (isotropic light)
    _ρ_1,_τ_1 = layer_ρ_τ(lha, bio, lwc, 1/MESOPHYLL_N, FT(90); N = N);
    _ρ_2,_τ_2 = layer_ρ_τ(lha, bio, lwc, 1 - 1/MESOPHYLL_N, FT(90); N = N);

    @. _denom = 1 - _ρ_1 * _ρ_2;
    @. _τ_leaf = _τ_s * _τ_2 / _denom;
    @. _ρ_leaf = _ρ_s + _τ_s * _ρ_2 * _τ_1 / _denom;

    # You can double check the conservation like this
    # _α_s = 1 .- _ρ_s .- _τ_s;
    # _α_1 = 1 .- _ρ_1 .- _τ_1;
    # _α_2 = 1 .- _ρ_2 .- _τ_2;
    # _α_leaf = _α_s .+ _τ_s .* _ρ_2 .* _α_1 ./ _denom .+ _τ_s .* _α_2 ./ _denom;

    return _ρ_leaf, _τ_leaf
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute the leaf level SIF matrices (before any reabsorption)
#
#######################################################################################################################################################################################################
"""

    leaf_raw_sif_matrices(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10) where {FT}

Return the leaf level SIF matrices (before any reabsorption), given
- `lha` leaf hyperspectral absorption coefficients
- `wls` hyperspectral wavelength set
- `bio` leaf hyperspectral biophysics
- `lwc` leaf water content
- `θ` incident angle of the incoming radiation
- `N` number of sublayers of the each sublayer (default: 10)

"""
function leaf_raw_sif_matrices(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10) where {FT}
    (; Φ_PS) = lha;
    (; IΛ_SIF, IΛ_SIFE) = wls;
    (; MESOPHYLL_N) = bio;

    # compute the reflectance and transmittance of the leaf layers (isotropic light)
    _ρ_s,_τ_s = layer_ρ_τ(lha, bio, lwc, 1/MESOPHYLL_N, θ; N = N);
    _ρ_1,_τ_1 = layer_ρ_τ(lha, bio, lwc, 1/MESOPHYLL_N, FT(90); N = N);
    _ρ_2,_τ_2 = layer_ρ_τ(lha, bio, lwc, 1 - 1/MESOPHYLL_N, FT(90); N = N);
    _denom = 1 .- _ρ_1 .* _ρ_2;

    # compute the mat_b_chl and mat_f_chl of the up and lower layers of a 2-layer leaf
    # _sife_s1 = layer_raw_sife(lha, wls, bio, lwc, 1/MESOPHYLL_N, θ; N = N);             # upper layer with light from upper side
    # _sife_l1 = layer_raw_sife(lha, wls, bio, lwc, 1/MESOPHYLL_N, FT(90); N = N);        # upper layer with light from lower side
    # _sife_l2 = layer_raw_sife(lha, wls, bio, lwc, 1 - 1/MESOPHYLL_N, FT(90); N = N);    # lower layer with light from upper side
    _sife_s1 = (1 .- _ρ_s .- _τ_s)[IΛ_SIFE];
    _sife_l1 = (1 .- _ρ_1 .- _τ_1)[IΛ_SIFE];
    _sife_l2 = (1 .- _ρ_2 .- _τ_2)[IΛ_SIFE];

    # compute the effective mat_b_chl and mat_f_chl of the leaf layers after counting for the SIF excitation wavelength
    # if you sum up _sife_1 .+ _sife_2, it should be equal to (1 - ρ_leaf - τ_leaf)[IΛ_SIFE]
    _sife_1 = _sife_s1 .+ _sife_l1 .* _τ_s[IΛ_SIFE] .* _ρ_2[IΛ_SIFE] ./ _denom[IΛ_SIFE];
    _sife_2 = _sife_l2 .* _τ_s[IΛ_SIFE] ./ _denom[IΛ_SIFE];

    # now compute the effective SIF emissions at the two layers (forward and backward)
    _mat_b_1 = _sife_1 ./ 2 * Φ_PS[IΛ_SIF]';
    _mat_f_1 = _sife_1 ./ 2 * Φ_PS[IΛ_SIF]';
    _mat_b_2 = _sife_2 ./ 2 * Φ_PS[IΛ_SIF]';
    _mat_f_2 = _sife_2 ./ 2 * Φ_PS[IΛ_SIF]';

    # here use the reflectance of a layer without absorption
    _ρ_zero,_τ_zero = layer_ρ_τ(lha, FT(90));
    _denom = 1 .- _ρ_zero .^ 2;

    # now consider the reflectance of SIF between the two layers to compute the matrices for the entire leaf (forward and backward)
    _mat_b = _mat_b_1 .+ _mat_f_1 .* _ρ_zero[IΛ_SIF]' .* _τ_zero[IΛ_SIF]' ./ _denom[IΛ_SIF]' .+ _mat_b_2 .* _τ_zero[IΛ_SIF]' ./ _denom[IΛ_SIF]';
    _mat_f = _mat_f_1 .* _τ_zero[IΛ_SIF]' ./ _denom[IΛ_SIF]' .+ _mat_f_2 .+ _mat_b_2 .* _ρ_zero[IΛ_SIF]' .* _τ_zero[IΛ_SIF]' ./ _denom[IΛ_SIF]';

    # the matrices should be conserved
    # @show sum(_sife_1 .+ _sife_2);
    # @show sum(_mat_b_1 .+ _mat_f_1 .+ _mat_b_2 .+ _mat_f_2);
    # @show sum(_mat_b .+ _mat_f)

    return _mat_b, _mat_f
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute the leaf level SIF matrices (after leaf reabsorption)
#
#######################################################################################################################################################################################################
"""

    leaf_sif_matrices(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10) where {FT}

Return the leaf level SIF matrices (after leaf reabsorption), given
- `lha` leaf hyperspectral absorption coefficients
- `wls` hyperspectral wavelength set
- `bio` leaf hyperspectral biophysics
- `lwc` leaf water content
- `θ` incident angle of the incoming radiation
- `N` number of sublayers of the each sublayer (default: 10)

"""
function leaf_sif_matrices(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10) where {FT}
    (; NR, Φ_PS) = lha;
    (; IΛ_SIF, IΛ_SIFE) = wls;
    (; MESOPHYLL_N) = bio;

    # reflectance and transmittance at the interface
    _τ_sf_sife = interface_isotropic_τ.(FT(1), NR[IΛ_SIFE], θ);
    _τ_21_sife = interface_isotropic_τ.(NR[IΛ_SIFE], FT(1), FT(90));
    _ρ_21_sife = 1 .- _τ_21_sife;
    _τ_12_sife = interface_isotropic_τ.(FT(1), NR[IΛ_SIFE], FT(90));

    # reflectance and transmittance through the layer
    _,_τ_0 = layer_ρ_τ(lha, bio, lwc, 1/MESOPHYLL_N, θ; N = N);
    _ρ_1,_τ_1 = layer_ρ_τ(lha, bio, lwc, 1/MESOPHYLL_N, FT(90); N = N);
    _ρ_2,_τ_2 = layer_ρ_τ(lha, bio, lwc, 1 - 1/MESOPHYLL_N, FT(90); N = N);

    #
    # PART 1: compute the matrices for the upper layer
    #

    # 1. compute the transmittance and absorption within a sublayer
    _t_sub = sublayer_τ(lha, bio, lwc, 1/MESOPHYLL_N, N);
    _t_sub_sife = _t_sub[IΛ_SIFE];
    _a_sub_sife = 1 .- _t_sub_sife;
    _t_all_sife = _t_sub_sife .^ N;
    _t_sub_sif  = _t_sub[IΛ_SIF];

    # 2. compute SIF emissions at the SIFE absorption site and rescale it based on the SIF absorption
    _rad_i = deepcopy(_τ_sf_sife);
    _mat_b_1 = zeros(FT, length(IΛ_SIFE), length(IΛ_SIF));
    _mat_f_1 = zeros(FT, length(IΛ_SIFE), length(IΛ_SIF));
    # now the radiation goes from up to down
    # the SIF up will be downscaled by _t_sub ^ (_i - 0.5) times
    # the SIF down will be downscaled by _t_sub ^ (N - _i + 0.5) times
    for _i in 1:N
        _mat_b_1 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^ (_i -0.5))';
        _mat_f_1 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^(N - _i + 0.5))';
        _rad_i .*= _t_sub_sife;
    end;
    # then the radiation goes from down to up
    _rad_i .*= _ρ_21_sife;
    for _i in N:-1:1
        _mat_b_1 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^ (_i -0.5))';
        _mat_f_1 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^(N - _i + 0.5))';
        _rad_i .*= _t_sub_sife;
    end;
    # then rescale it by 1 / (1 - ρ_21 * t_all * ρ_21 * t_all) for the SIFE
    _denom = 1 .- _t_all_sife .* _ρ_21_sife .* _t_all_sife .* _ρ_21_sife;
    _mat_b_1 ./= _denom;
    _mat_f_1 ./= _denom;

    #
    # Part 2: compute the matrices for the lower layer
    #

    # 1. compute the transmittance and absorption within a sublayer
    _t_sub = sublayer_τ(lha, bio, lwc, 1 - 1/MESOPHYLL_N, N);
    _t_sub_sife = _t_sub[IΛ_SIFE];
    _a_sub_sife = 1 .- _t_sub_sife;
    _t_all_sife = _t_sub_sife .^ N;
    _t_sub_sif  = _t_sub[IΛ_SIF];

    # 2. compute SIF emissions at the SIFE absorption site and rescale it based on the SIF absorption
    _rad_i = _τ_0[IΛ_SIFE] .* _τ_12_sife;
    _mat_b_2 = zeros(FT, length(IΛ_SIFE), length(IΛ_SIF));
    _mat_f_2 = zeros(FT, length(IΛ_SIFE), length(IΛ_SIF));
    # now the radiation goes from up to down
    for _i in 1:N
        _mat_b_2 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^ (_i -0.5))';
        _mat_f_2 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^(N - _i + 0.5))';
        _rad_i .*= _t_sub_sife;
    end;
    # then the radiation goes from down to up
    _rad_i .*= _ρ_21_sife;
    for _i in N:-1:1
        _mat_b_2 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^ (_i -0.5))';
        _mat_f_2 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^(N - _i + 0.5))';
        _rad_i .*= _t_sub_sife;
    end;
    # then rescale it by 1 / (1 - ρ_21 * t_all * ρ_21 * t_all) for the SIFE
    _denom = 1 .- _t_all_sife .* _ρ_21_sife .* _t_all_sife .* _ρ_21_sife;
    _mat_b_2 ./= _denom;
    _mat_f_2 ./= _denom;

    #
    # Part 3: compute the matrices for the entire leaf
    #

    # 1. here use the reflectance of a layer with absorption for SIF
    _ρ_1_sif = _ρ_1[IΛ_SIF];
    _τ_1_sif = _τ_1[IΛ_SIF];
    _ρ_2_sif = _ρ_2[IΛ_SIF];
    _τ_2_sif = _τ_2[IΛ_SIF];
    _denom = 1 .- _ρ_1_sif .* _ρ_2_sif;

    # now consider the reflectance and transmittance of SIF between the two layers to compute the matrices for the entire leaf (forward and backward)
    _mat_b = _mat_b_1 .+ _mat_f_1 .* _ρ_2_sif' .* _τ_1_sif' ./ _denom' .+ _mat_b_2 .* _τ_1_sif' ./ _denom';
    _mat_f = _mat_f_1 .* _τ_2_sif' ./ _denom' .+ _mat_f_2 .+ _mat_b_2 .* _ρ_1_sif' .* _τ_2_sif' ./ _denom';

    return _mat_b, _mat_f
end;
