#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2025-Feb-09: add method for Dualspect as an extended Fluspect with inter-plate scattering
# Bug fix
#     2025-May-20: fix an bug related to the temporary mat_b and mat_f for the two layers
#
#######################################################################################################################################################################################################
leaf_sif_matrices!(config::SPACConfiguration{FT}, bio::LeafBio{FT}, mtd::SIFMatrixDualspectMethod) where {FT} = (
    (; IΛ_SIF, IΛ_SIFE) = config.SPECTRA;
    (; ρ_interface_θ, τ_interface_θ, ρ_interface_12, τ_interface_12, ρ_interface_21, τ_interface_21, f_sife) = bio.auxil;
    (; ρ_interface_12_eff, τ_interface_12_eff, ρ_interface_21_eff, τ_interface_21_eff) = bio.auxil;
    (; ρ_layer_θ, τ_layer_θ, ρ_layer_1, τ_layer_1, ρ_layer_2, τ_layer_2) = bio.auxil;

    # 1. Get the matrix for the first layer (1 plate) using incoming radiation angle
    mat_b_θ = similar(bio.auxil.mat_b);
    mat_f_θ = similar(bio.auxil.mat_f);
    kubelka_munk_sif_matrices!(config, ρ_layer_θ, τ_layer_θ, ρ_interface_θ, τ_interface_θ, ρ_interface_21, τ_interface_21, f_sife, mtd.N, mat_b_θ, mat_f_θ);

    # 2. Get the matrix for the first layer (1 plate) using diffuse radiation
    mat_b_1 = similar(bio.auxil.mat_b);
    mat_f_1 = similar(bio.auxil.mat_f);
    kubelka_munk_sif_matrices!(config, ρ_layer_1, τ_layer_1, ρ_interface_12, τ_interface_12, ρ_interface_21, τ_interface_21, f_sife, mtd.N, mat_b_1, mat_f_1);

    # 3. Get the matrix for the second layer (N-1 plates) using diffuse radiation
    mat_b_2 = similar(bio.auxil.mat_b);
    mat_f_2 = similar(bio.auxil.mat_f);
    kubelka_munk_sif_matrices!(config, ρ_layer_2, τ_layer_2, ρ_interface_12_eff, τ_interface_12_eff, ρ_interface_21_eff, τ_interface_21_eff, f_sife, mtd.N, mat_b_2, mat_f_2);

    # 4. Now that we need to consider 3 cases for non-scattered SIF:
    #     - case 1: incoming radiation directly to the first layer
    #     - case 2: incoming radiation that is transmitted and then scattered to the second layer
    #     - case 3: incoming radiation that is transmitted and then scattered to the second layer, and then reflected back to the first layer
    # The radiation profiles are
    #     - case 1: 1
    #     - case 2: τ_1 / (1 - ρ_2 * ρ_1)
    #     - case 3: τ_1 * ρ_2 / (1 - ρ_2 * ρ_1)
    # Therefore,
    #     the  forward SIF matrix of the  first layer is sum of case 1  forward SIF and case 3 backward SIF
    #     the backward SIF matrix of the  first layer is sum of case 1 backward SIF and case 2  forward SIF
    #     the forward/backward SIF matrices of the second layer are scaled based on the radiation
    rad_2 = @. τ_layer_1 / (1 - ρ_layer_2 * ρ_layer_1);
    rad_3 = @. τ_layer_1 * ρ_layer_2 / (1 - ρ_layer_2 * ρ_layer_1);
    rad_2_sife = rad_2[IΛ_SIFE];
    rad_3_sife = rad_3[IΛ_SIFE];
    _mat_b_1 = similar(bio.auxil.mat_b);
    _mat_f_1 = similar(bio.auxil.mat_f);
    _mat_b_2 = similar(bio.auxil.mat_b);
    _mat_f_2 = similar(bio.auxil.mat_f);
    _mat_b_1 .= mat_b_θ .+ mat_f_1 .* rad_3_sife';
    _mat_f_1 .= mat_f_θ .+ mat_b_1 .* rad_3_sife';
    _mat_b_2 .= mat_b_2 .* rad_2_sife';
    _mat_f_2 .= mat_f_2 .* rad_2_sife';

    # 5. Now, consider the SIF scattering between the 2 layers.
    # The Total backward SIF is the sum of
    #     - the backward SIF of the  first layer
    #     - the  forward SIF of the  first layer that got scattered and transmitted
    #     - the backward SIF of the second layer that got transmitted and scattered
    # The Total forward SIF is the sum of
    #     - the  forward SIF of the  first layer that got scattered and transmitted
    #     - the backward SIF of the second layer that got scattered and transmitted
    #     - the  forward SIF of the second layer
    ρ_1_sif = ρ_layer_1[IΛ_SIF];
    τ_1_sif = τ_layer_1[IΛ_SIF];
    ρ_2_sif = ρ_layer_2[IΛ_SIF];
    τ_2_sif = τ_layer_2[IΛ_SIF];
    denom = 1 .- ρ_2_sif .* ρ_1_sif;
    bio.auxil.mat_b .= _mat_b_1 .+ _mat_f_1 .* ρ_2_sif .* τ_1_sif ./ denom .+ _mat_b_2 .* τ_1_sif ./ denom;
    bio.auxil.mat_f .= _mat_f_1 .* τ_2_sif ./ denom .+ _mat_b_2 .* ρ_1_sif .* τ_2_sif ./ denom .+ _mat_f_2;

    # compute the mean and mean diff of mat_b and mat_f
    bio.auxil.mat_mean .= (bio.auxil.mat_b .+ bio.auxil.mat_f) ./ 2;
    bio.auxil.mat_diff .= (bio.auxil.mat_b .- bio.auxil.mat_f) ./ 2;

    return nothing
);
