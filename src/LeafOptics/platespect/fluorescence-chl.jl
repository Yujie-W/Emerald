"""

    layer_1_sif_vec_chl!(τ_i_θ::FT, τ_i_12::FT, τ_i_21::FT, k_e::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, f_sife::FT, vec_b::SubArray, vec_f::SubArray) where {FT}

Update the SIF conversion matrix of the first layer without reabsorption, given
- `τ_i_θ` transmittance of the incoming radiation at the air-water interface
- `τ_i_12` transmittance of the isotropic radiation at the air-water interface
- `τ_i_21` transmittance of the isotropic radiation at the water-air interface
- `k_e` extinction coefficient within a sublayer of layer 1
- `τ_θ` transmittance of the incoming radiation across the leaf layer
- `ρ_1` reflectance of the first layer (1)
- `ρ_2` reflectance of the second layer (n-1)
- `f_sife` SIF excitation scaling factor (f_cab + f_car * ϕ_car)
- `vec_b` SIF vector backwards (SubArray of a Matrix)
- `vec_f` SIF vector forwards (SubArray of a Matrix)

"""
function layer_1_sif_vec_chl!(τ_i_θ::FT, τ_i_12::FT, τ_i_21::FT, k_e::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, f_sife::FT, vec_b::SubArray, vec_f::SubArray) where {FT}
    # parameters required for the calculation that can be derived from the input parameters
    ρ_i_21 = 1 - τ_i_21;
    τ_all = exp(-k_e);
    denom = 1 - τ_all * ρ_i_21 * τ_all * ρ_i_21;

    # computing e_up and e_down
    e_1ꜜ = τ_i_θ / denom;
    e_1ꜛ = e_1ꜜ * τ_all * ρ_i_21;
    e_2ꜛ = τ_θ * ρ_2 / (1 - ρ_1 * ρ_2) * τ_i_12 / denom;
    e_2ꜜ = e_2ꜛ * τ_all * ρ_i_21;
    eꜜ = e_1ꜜ + e_2ꜜ;
    eꜛ = e_1ꜛ + e_2ꜛ;

    # compute the SIF emission spectrum with integrated equation
    @. vec_f = (eꜜ + eꜛ) * (1 - τ_all) * f_sife / 2;
    @. vec_b = (eꜜ + eꜛ) * (1 - τ_all) * f_sife / 2;

    return nothing
end;


"""

    layer_2_sif_vec_chl!(k_e::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, τ_2::FT, f_sife::FT, vec_b::SubArray, vec_f::SubArray) where {FT}

Update the SIF conversion matrix of the n-1 layer without reabsorption, given
- `k_e` extinction coefficient within a sublayer of layer 2
- `τ_θ` transmittance of the incoming radiation across the leaf layer
- `ρ_1` reflectance of the first layer (1)
- `ρ_2` reflectance of the second layer (n-1)
- `τ_2` transmittance of the second layer (n-1)
- `f_sife` SIF excitation scaling factor (f_cab + f_car * ϕ_car)
- `vec_b` SIF vector backwards (SubArray of a Matrix)
- `vec_f` SIF vector forwards (SubArray of a Matrix)

"""
function layer_2_sif_vec_chl!(k_e::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, τ_2::FT, f_sife::FT, vec_b::SubArray, vec_f::SubArray) where {FT}
    # parameters required for the calculation
    τ_all = exp(-k_e);

    # 1. here we consider the n-1 layers as one single layer, and the SIF transmission within this effective layer is same as the computed τ_sub_2
    #    then we need to rescale the interface ρ and τ for the effective layer so that the computed layer level ρ and τ are same as computed
    ρ_i_12 = effective_ρ_12(ρ_2, τ_2, τ_all);
    ρ_i_21 = effective_ρ_21(ρ_2, τ_2, τ_all);
    τ_i_12 = 1 - ρ_i_12;
    denom = 1 - τ_all * ρ_i_21 * τ_all * ρ_i_21;

    # computing e_up and e_down
    eꜜ = τ_θ / (1 - ρ_1 * ρ_2) * τ_i_12 / denom;
    eꜛ = eꜜ * τ_all * ρ_i_21;

    # compute the SIF emission spectrum  with integrated equation
    @. vec_f = (eꜜ + eꜛ) * (1 - τ_all) * f_sife / 2;
    @. vec_b = (eꜜ + eꜛ) * (1 - τ_all) * f_sife / 2;

    return nothing
end;


"""

    leaf_sif_matrices_chl!(config::SPACConfig{FT}, bio::LeafBio{FT}, cache::SPACCache{FT}) where {FT}

Update the SIF conversion matrix of the leaf without reabsorption, given
- `config` SPAC configuration
- `bio` leaf biophysics

"""
function leaf_sif_matrices_chl! end;

leaf_sif_matrices_chl!(config::SPACConfig{FT}, bio::LeafBio{FT}, cache::SPACCache{FT}) where {FT} = leaf_sif_matrices_chl!(config, bio, cache, config.METHODS.FLUORESCENCE_SPECTRA_METHOD);

leaf_sif_matrices_chl!(config::SPACConfig{FT}, bio::LeafBio{FT}, cache::SPACCache{FT}, ::PlatespectFluorescenceSpectra) where {FT} = (
    (; SPECTRA) = config.CONSTANTS;
    (; IΛ_SIF, IΛ_SIFE, ΔΛ_SIF, Λ_SIF, Λ_SIFE, Φ_PS) = SPECTRA;

    # update the SIF emission vector per excitation wavelength
    ϕ           = bio.auxil._ϕ_sif;
    factor      = cache.cache_sif_1;

    for i in eachindex(IΛ_SIFE)
        ii = IΛ_SIFE[i];

        # read the SIF emission spectrum
        ϕ .= view(Φ_PS, IΛ_SIF);

        # tune SIF emission PDF based on the SIF excitation wavelength
        expsife = exp(Λ_SIFE[ii] / 10);
        @. factor = 1 / (1 + exp(-Λ_SIF / 10) * expsife);
        ϕ .*= factor;

        # rescale ϕ
        ϕ ./= ΔΛ_SIF' * ϕ;

        # read in the values from the auxiliary variables
        vec_b_1     = view(bio.auxil.mat_b_1_chl, :, i);
        vec_f_1     = view(bio.auxil.mat_f_1_chl, :, i);
        vec_b_2     = view(bio.auxil.mat_b_2_chl, :, i);
        vec_f_2     = view(bio.auxil.mat_f_2_chl, :, i);

        vec_b       = view(bio.auxil.mat_b_chl, :, i);
        vec_f       = view(bio.auxil.mat_f_chl, :, i);

        τ_i_θ       = bio.auxil.τ_interface_θ[ii];      # the transmittance of the incoming radiation at the air-water interface
        τ_i_12      = bio.auxil.τ_interface_12[ii];     # the transmittance of the isotropic radiation at the air-water interface
        τ_i_21      = bio.auxil.τ_interface_21[ii];     # the transmittance of the isotropic radiation at the water-air interface
        k_all_1     = bio.auxil.k_all_1[ii];            # the extinction coefficient within a sublayer of layer 1
        k_all_2     = bio.auxil.k_all_2[ii];            # the extinction coefficient within a sublayer of layer 2 (n-1)
        τ_l_θ       = bio.auxil.τ_layer_θ[ii];          # the transmittance of the incoming radiation across the leaf layer 1
        ρ_l_1       = bio.auxil.ρ_layer_1[ii];          # the reflectance of isotropic radiation across layer 1
        ρ_l_2       = bio.auxil.ρ_layer_2[ii];          # the reflectance of isotropic radiation across layer 2 (n-1)
        τ_l_2       = bio.auxil.τ_layer_2[ii];          # the transmittance of isotropic radiation across layer 2 (n-1)
        f_sife      = bio.auxil.f_sife[ii];

        # update the SIF conversion matrix of the two layers (SIF that reachs the internal the water-air interface)
        layer_1_sif_vec_chl!(τ_i_θ, τ_i_12, τ_i_21, k_all_1, τ_l_θ, ρ_l_1, ρ_l_2, f_sife, vec_b_1, vec_f_1);
        layer_2_sif_vec_chl!(k_all_2, τ_l_θ, ρ_l_1, ρ_l_2, τ_l_2, f_sife, vec_b_2, vec_f_2);

        # compute the SIF emission vector backward and forward
        vec_b .= vec_b_1 .+ vec_b_2;
        vec_f .= vec_f_1 .+ vec_f_2;

        # scale the matrices based on the Φ_PS*
        vec_b .*= ϕ;
        vec_f .*= ϕ;
    end;

    # compute the mean and mean diff of mat_b_chl and mat_f_chl
    bio.auxil.matꜛ_chl .= (bio.auxil.mat_b_chl .+ bio.auxil.mat_f_chl) ./ 2;
    bio.auxil.matꜜ_chl .= (bio.auxil.mat_b_chl .- bio.auxil.mat_f_chl) ./ 2;

    return nothing
);
