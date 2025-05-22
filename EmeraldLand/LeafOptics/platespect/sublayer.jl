#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2025-05-22: adding new function for PlateSpect method (may be easier to understand)
#
#######################################################################################################################################################################################################
function radiation_triggered_sif_vec!(τ_i_θ::FT, τ_i_21::FT, k_e::FT, f_sife::FT, k_f::Vector{FT}, vec_b::Vector{FT}, vec_f::Vector{FT}) where {FT}
    # parameters required for the calculation that can be derived from the input parameters
    ρ_i_21 = 1 - τ_i_21;
    τ_all = exp(-k_e);
    denom = 1 - τ_all * ρ_i_21 * τ_all * ρ_i_21;

    # computing e_up and e_down
    eꜜ = τ_i_θ / denom;
    eꜛ = eꜜ * τ_all * ρ_i_21;

    # vec b and vec f are the SIF emission spectrum for the backward and forward directions respectively
    normꜜ = k_e * f_sife / 2 * eꜜ;
    normꜛ = k_e * f_sife / 2 * eꜛ;

    @. vec_f = normꜜ * int_exp_kef⁻(k_e,k_f) + normꜛ * int_exp_kef⁺(k_e, k_f);
    @. vec_b = normꜛ * int_exp_kef⁻(k_e,k_f) + normꜜ * int_exp_kef⁺(k_e, k_f);

    return nothing
end;


function leaf_sif_matrices_new!(config::SPACConfiguration{FT}, bio::LeafBio{FT}, ::SIFMatrixPlatespectMethod) where {FT}
    (; Φ_SIF_CUTOFF, Φ_SIF_RESCALE) = config;
    (; IΛ_SIF, IΛ_SIFE, ΔΛ_SIF, Λ_SIF, Λ_SIFE, Φ_PS) = config.SPECTRA;

    # alias from the auxiliary variables
    ρ_i_21     = bio.auxil.ρ_interface_21;
    τ_i_θ      = bio.auxil.τ_interface_θ;
    τ_i_12     = bio.auxil.τ_interface_12;
    τ_i_21     = bio.auxil.τ_interface_21;
    ρ_i_21_eff = bio.auxil.ρ_interface_21_eff;
    τ_i_12_eff = bio.auxil.τ_interface_12_eff;
    τ_i_21_eff = bio.auxil.τ_interface_21_eff;
    k_all_1    = bio.auxil.k_all_1;
    k_all_2    = bio.auxil.k_all_2;
    τ_all_1    = bio.auxil.τ_all_1;
    τ_all_2    = bio.auxil.τ_all_2;
    ρ_layer_1  = bio.auxil.ρ_layer_1;
    ρ_layer_2  = bio.auxil.ρ_layer_2;
    τ_layer_1  = bio.auxil.τ_layer_1;
    τ_layer_2  = bio.auxil.τ_layer_2;

    # for fluorescence
    k_f_1 = k_all_1[IΛ_SIF];
    k_f_2 = k_all_2[IΛ_SIF];
    τ_f_1 = τ_all_1[IΛ_SIF];
    τ_f_2 = τ_all_2[IΛ_SIF];
    ρ_i_21_f_1 = ρ_i_21[IΛ_SIF];
    ρ_i_21_f_2 = ρ_i_21_eff[IΛ_SIF];
    ρ_1_f = ρ_layer_1[IΛ_SIF];
    ρ_2_f = ρ_layer_2[IΛ_SIF];
    τ_1_f = τ_layer_1[IΛ_SIF];
    τ_2_f = τ_layer_2[IΛ_SIF];

    f_b_θ = similar(k_f_1);
    f_f_θ = similar(k_f_1);
    f_b_1 = similar(k_f_1);
    f_f_1 = similar(k_f_1);
    f_b_2 = similar(k_f_1);
    f_f_2 = similar(k_f_1);

    ϕ = bio.auxil._ϕ_sif;

    for i in eachindex(IΛ_SIFE)
        i_e = IΛ_SIFE[i];

        # compute the radiation triggered without scaling or scattering
        k_e_1 = k_all_1[i_e];
        k_e_2 = k_all_2[i_e];
        f_sife = bio.auxil.f_sife[i_e];
        τ_i_12_e_θ = τ_i_θ[i_e];
        τ_i_12_e_1 = τ_i_12[i_e];
        τ_i_21_e_1 = τ_i_21[i_e];
        τ_i_12_e_2 = τ_i_12_eff[i_e];
        τ_i_21_e_2 = τ_i_21_eff[i_e];
        radiation_triggered_sif_vec!(τ_i_12_e_θ, τ_i_21_e_1, k_e_1, f_sife, k_f_1, f_b_θ, f_f_θ);
        radiation_triggered_sif_vec!(τ_i_12_e_1, τ_i_21_e_1, k_e_1, f_sife, k_f_1, f_b_1, f_f_1);
        radiation_triggered_sif_vec!(τ_i_12_e_2, τ_i_21_e_2, k_e_2, f_sife, k_f_2, f_b_2, f_f_2);

        # scale the radiation triggered SIF by the radiation intensity
        ρ_1_e = bio.auxil.ρ_layer_1[i_e];
        ρ_2_e = bio.auxil.ρ_layer_2[i_e];
        τ_1_e = bio.auxil.τ_layer_1[i_e];
        rad_θ = 1;
        rad_1 = τ_1_e / (1 - ρ_1_e * ρ_2_e) * ρ_2_e;
        rad_2 = τ_1_e / (1 - ρ_1_e * ρ_2_e);
        @. f_b_θ *= rad_θ;
        @. f_f_θ *= rad_θ;
        @. f_b_1 *= rad_1;
        @. f_f_1 *= rad_1;
        @. f_b_2 *= rad_2;
        @. f_f_2 *= rad_2;

        # add the triggered SIF up that reached the interface (not rescatter here yet)
        f_i_b_1 = @. f_b_θ + f_f_1;
        f_i_f_1 = @. f_f_θ + f_b_1;
        f_i_b_2 = @. f_b_2;
        f_i_f_2 = @. f_f_2;

        # compute the SIF out of the interface
        f_out_b_1 = @. layer_sif_out(f_i_f_1, f_i_b_1, ρ_i_21_f_1, τ_f_1);
        f_out_f_1 = @. layer_sif_out(f_i_b_1, f_i_f_1, ρ_i_21_f_1, τ_f_1);
        f_out_b_2 = @. layer_sif_out(f_i_f_2, f_i_b_2, ρ_i_21_f_2, τ_f_2);
        f_out_f_2 = @. layer_sif_out(f_i_b_2, f_i_f_2, ρ_i_21_f_2, τ_f_2);

        # compute the SIF out of the leaf
        f_b = @. f_out_b_1 +
                 f_out_f_1 * ρ_2_f * τ_1_f / (1 - ρ_1_f * ρ_2_f) +
                 f_out_b_2 * τ_1_f / (1 - ρ_1_f * ρ_2_f);
        f_f = @. f_out_f_1 * τ_2_f / (1 - ρ_1_f * ρ_2_f) +
                 f_out_b_2 * ρ_1_f * τ_2_f / (1 - ρ_1_f * ρ_2_f) +
                 f_out_f_2;

        # scale the fluorescence by the PDF
        ϕ .= view(Φ_PS, IΛ_SIF);
        #  tune SIF emission PDF based on the SIF excitation wavelength
        #     0: not cut off
        #     1: sharp cut off
        #     2: sigmoid cut off used in SCOPE
        #     x: add more functions if you wish
        if Φ_SIF_CUTOFF == 1
            if i_e > IΛ_SIF[1]
                ϕ[1:i_e-IΛ_SIF[1]] .= 0;
            end;
        elseif Φ_SIF_CUTOFF == 2
            factor = 1 ./ (1 .+ exp.(-Λ_SIF ./ 10) .* exp(Λ_SIFE[i] / 10));
            ϕ .*= factor;
        end;
         # rescale ϕ if Φ_SIF_RESCALE is true
        if Φ_SIF_RESCALE && Φ_SIF_CUTOFF > 0
            ϕ ./= ΔΛ_SIF' * ϕ;
        end;
        @. f_b *= ϕ;
        @. f_f *= ϕ;

        # update the matrices
        @. bio.auxil.mat_b[:,i_e] .= f_b;
        @. bio.auxil.mat_f[:,i_e] .= f_f;
    end;

    # compute the mean and mean diff of mat_b and mat_f
    bio.auxil.mat_mean .= (bio.auxil.mat_b .+ bio.auxil.mat_f) ./ 2;
    bio.auxil.mat_diff .= (bio.auxil.mat_b .- bio.auxil.mat_f) ./ 2;

    return nothing
end;
