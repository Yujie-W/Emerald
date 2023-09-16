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


#=
#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute the leaf level SIF matrices (before any reabsorption)
#     2023-Sep-14: use layer_ρ_τ_direct and layer_ρ_τ_diffuse to compute the reflectance and transmittance of the leaf layers
#     2023-Sep-14: use only cab and car absorption for the SIF excitation
#
#######################################################################################################################################################################################################
"""

    leaf_raw_sif_matrices(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10, ϕ_car::FT = FT(0)) where {FT}

Return the leaf level SIF matrices (before any reabsorption), given
- `lha` leaf hyperspectral absorption coefficients
- `wls` hyperspectral wavelength set
- `bio` leaf hyperspectral biophysics
- `lwc` leaf water content
- `θ` incident angle of the incoming radiation
- `N` number of sublayers of each layer
- `ϕ_car` carotenoid contribution to chlorophyll fluorescence (default: 0)

"""
function leaf_raw_sif_matrices(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10, ϕ_car::FT = FT(0)) where {FT}
    (; Φ_PS) = lha;
    (; IΛ_SIF, IΛ_SIFE) = wls;
    (; MESOPHYLL_N) = bio;

    # compute the reflectance and transmittance of the leaf layers (isotropic light)
    # _ρ_s,_τ_s = layer_ρ_τ(lha, bio, lwc, 1/MESOPHYLL_N, θ; N = N);
    # _ρ_1,_τ_1 = layer_ρ_τ(lha, bio, lwc, 1/MESOPHYLL_N, FT(90); N = N);
    # _ρ_2,_τ_2 = layer_ρ_τ(lha, bio, lwc, 1 - 1/MESOPHYLL_N, FT(90); N = N);
    _ρ_s,_τ_s = layer_ρ_τ_direct(lha, bio, lwc, θ; N= N);
    _ρ_1,_τ_1,_ρ_2,_τ_2 = layer_ρ_τ_diffuse(lha, bio, lwc; N = N);
    _denom = 1 .- _ρ_1 .* _ρ_2;

    # compute the absorption by cab and car
    _,_f_cab,_f_car = sublayer_τ(lha, bio, lwc, 1/MESOPHYLL_N, N);
    _ϕ_sife = (_f_cab .+ _f_car .* ϕ_car)[IΛ_SIFE];

    # compute the mat_b_chl and mat_f_chl of the up and lower layers of a 2-layer leaf
    # _sife_s1 = layer_raw_sife(lha, wls, bio, lwc, 1/MESOPHYLL_N, θ; N = N);             # upper layer with light from upper side
    # _sife_l1 = layer_raw_sife(lha, wls, bio, lwc, 1/MESOPHYLL_N, FT(90); N = N);        # upper layer with light from lower side
    # _sife_l2 = layer_raw_sife(lha, wls, bio, lwc, 1 - 1/MESOPHYLL_N, FT(90); N = N);    # lower layer with light from upper side
    _sife_s1 = (1 .- _ρ_s .- _τ_s)[IΛ_SIFE] .* _ϕ_sife;
    _sife_l1 = (1 .- _ρ_1 .- _τ_1)[IΛ_SIFE] .* _ϕ_sife;
    _sife_l2 = (1 .- _ρ_2 .- _τ_2)[IΛ_SIFE] .* _ϕ_sife;

    # compute the effective mat_b_chl and mat_f_chl of the leaf layers after counting for the SIF excitation wavelength
    # if you sum up _sife_1 .+ _sife_2, it should be equal to (1 - ρ_leaf - τ_leaf)[IΛ_SIFE]
    _sife_1 = _sife_s1 .+ _sife_l1 .* _τ_s[IΛ_SIFE] .* _ρ_2[IΛ_SIFE] ./ _denom[IΛ_SIFE];
    _sife_2 = _sife_l2 .* _τ_s[IΛ_SIFE] ./ _denom[IΛ_SIFE];

    # now compute the effective SIF emissions at the two layers (forward and backward)
    # TODO: here the assumption that SIF is evenly distributed in the lower layer is imperfect because the lower n-1 layers could have internal reflectance and transmittance
    #       a better way would be to modify the the leaf_sif_matrices function
    #       currently, this is enough to compute the total SIF emissions from the chlorophyll
    _mat_b_1 = _sife_1 ./ 2 * Φ_PS[IΛ_SIF]';
    _mat_f_1 = _sife_1 ./ 2 * Φ_PS[IΛ_SIF]';
    _mat_b_2 = _sife_2 ./ 2 * Φ_PS[IΛ_SIF]';
    _mat_f_2 = _sife_2 ./ 2 * Φ_PS[IΛ_SIF]';

    # here use the reflectance of a layer without absorption
    #_ρ_zero,_τ_zero = layer_ρ_τ(lha, FT(90));
    _ρ0_1,_τ0_1,_ρ0_2,_τ0_2 = layer_ρ_τ_diffuse(lha, bio);
    _denom = 1 .- _ρ0_1 .* _ρ0_2;

    # now consider the reflectance of SIF between the two layers to compute the matrices for the entire leaf (forward and backward)
    _mat_b = _mat_b_1 .+ _mat_f_1 .* _ρ0_2[IΛ_SIF]' .* _τ0_1[IΛ_SIF]' ./ _denom[IΛ_SIF]' .+ _mat_b_2 .* _τ0_1[IΛ_SIF]' ./ _denom[IΛ_SIF]';
    _mat_f = _mat_f_1 .* _τ0_2[IΛ_SIF]' ./ _denom[IΛ_SIF]' .+ _mat_f_2 .+ _mat_b_2 .* _ρ0_1[IΛ_SIF]' .* _τ0_2[IΛ_SIF]' ./ _denom[IΛ_SIF]';

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
#     2023-Sep-14: use layer_ρ_τ_direct and layer_ρ_τ_diffuse to compute the reflectance and transmittance of the leaf layers
#     2023-Sep-14: use only cab and car absorption for the SIF excitation
#
#######################################################################################################################################################################################################
"""

    leaf_sif_matrices(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10, ϕ_car::FT = FT(0)) where {FT}

Return the leaf level SIF matrices (after leaf reabsorption), given
- `lha` leaf hyperspectral absorption coefficients
- `wls` hyperspectral wavelength set
- `bio` leaf hyperspectral biophysics
- `lwc` leaf water content
- `θ` incident angle of the incoming radiation
- `N` number of sublayers of each layer
- `ϕ_car` carotenoid contribution to chlorophyll fluorescence (default: 0)

"""
function leaf_sif_matrices(lha::HyperspectralAbsorption{FT}, wls::WaveLengthSet{FT}, bio::HyperspectralLeafBiophysics{FT}, lwc::FT, θ::FT = FT(40); N::Int = 10, ϕ_car::FT = FT(0)) where {FT}
    (; NR, Φ_PS) = lha;
    (; IΛ_SIF, IΛ_SIFE) = wls;
    (; MESOPHYLL_N) = bio;

    # reflectance and transmittance at the interface
    _τ_sf_sife = interface_isotropic_τ.(FT(1), NR[IΛ_SIFE], θ);
    _τ_21_sife = interface_isotropic_τ.(NR[IΛ_SIFE], FT(1), FT(90));
    _ρ_21_sife = 1 .- _τ_21_sife;

    # reflectance and transmittance through the layer
    # _,_τ_0 = layer_ρ_τ(lha, bio, lwc, 1/MESOPHYLL_N, θ; N = N);
    # _ρ_1,_τ_1 = layer_ρ_τ(lha, bio, lwc, 1/MESOPHYLL_N, FT(90); N = N);
    # _ρ_2,_τ_2 = layer_ρ_τ(lha, bio, lwc, 1 - 1/MESOPHYLL_N, FT(90); N = N);
    _,_τ_0 = layer_ρ_τ_direct(lha, bio, lwc, θ; N= N);
    _ρ_1,_τ_1,_ρ_2,_τ_2 = layer_ρ_τ_diffuse(lha, bio, lwc; N = N);

    #
    # PART 1: compute the matrices for the upper layer
    #

    # 1. compute the transmittance and absorption within a sublayer
    _t_sub,_f_cab,_f_car = sublayer_τ(lha, bio, lwc, 1/MESOPHYLL_N, N);
    _t_sub_sife = _t_sub[IΛ_SIFE];
    _a_sub_sife = 1 .- _t_sub_sife;
    _t_all_sife = _t_sub_sife .^ N;
    _t_sub_sif  = _t_sub[IΛ_SIF];
    _ϕ_sife = (_f_cab .+ _f_car .* ϕ_car)[IΛ_SIFE];

    # 2. compute SIF emissions at the SIFE absorption site and rescale it based on the SIF absorption
    _rad_i = _τ_sf_sife .* _ϕ_sife;
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
    _t_sub,_,_ = sublayer_τ(lha, bio, lwc, 1 - 1/MESOPHYLL_N, N);
    _t_sub_sife = _t_sub[IΛ_SIFE];
    _a_sub_sife = 1 .- _t_sub_sife;
    _t_all_sife = _t_sub_sife .^ N;
    _t_sub_sif  = _t_sub[IΛ_SIF];

    # 2. compute the effective τ_21 from the effective reflectance and transmittance of the n-1 layers
    _τ_21_rescale = (_t_all_sife .* _τ_2[IΛ_SIFE] .+ 1 .- _ρ_2[IΛ_SIFE]) ./ (_t_all_sife .* _τ_2[IΛ_SIFE] .+ NR[IΛ_SIFE]);
    _τ_12_rescale = _τ_21_rescale .* NR[IΛ_SIFE];
    _ρ_21_rescale = 1 .- _τ_21_rescale;

    # 3. compute SIF emissions at the SIFE absorption site and rescale it based on the SIF absorption
    _rad_i = _τ_0[IΛ_SIFE] .* _ϕ_sife .* _τ_12_rescale;
    _mat_b_2 = zeros(FT, length(IΛ_SIFE), length(IΛ_SIF));
    _mat_f_2 = zeros(FT, length(IΛ_SIFE), length(IΛ_SIF));
    # now the radiation goes from up to down
    for _i in 1:N
        _mat_b_2 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^ (_i -0.5))';
        _mat_f_2 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^(N - _i + 0.5))';
        _rad_i .*= _t_sub_sife;
    end;
    # then the radiation goes from down to up
    _rad_i .*= _ρ_21_rescale;
    for _i in N:-1:1
        _mat_b_2 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^ (_i -0.5))';
        _mat_f_2 .+= (_rad_i ./ 2 .* _a_sub_sife) * Φ_PS[IΛ_SIF]' .* (_t_sub_sif .^(N - _i + 0.5))';
        _rad_i .*= _t_sub_sife;
    end;
    # then rescale it by 1 / (1 - ρ_21 * t_all * ρ_21 * t_all) for the SIFE
    _denom = 1 .- _t_all_sife .* _ρ_21_rescale .* _t_all_sife .* _ρ_21_rescale;
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
=#
