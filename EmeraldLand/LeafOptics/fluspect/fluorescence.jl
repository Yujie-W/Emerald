#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-21: add function to compute the matrices using doubling adding method
#     2024-Aug-16: refactor the function to use SIFMatrixFluspectMethod (N_DUB)
#     2025-Feb-09: generalize the double adding method to work for both Fluspect and Dualspect
#
#######################################################################################################################################################################################################
leaf_sif_matrices!(config::SPACConfiguration{FT}, bio::LeafBio{FT}, mtd::SIFMatrixFluspectMethod) where {FT} = (
    (; ρ_leaf, τ_leaf, ρ_interface_θ, τ_interface_θ, ρ_interface_21, τ_interface_21, f_sife, mat_b, mat_f) = bio.auxil;

    # update the mat_b and mat_f based on the doubling adding method (which has been generalized to work for both Fluspect and Dualspect)
    kubelka_munk_sif_matrices_new!(config, ρ_leaf, τ_leaf, ρ_interface_θ, τ_interface_θ, ρ_interface_21, τ_interface_21, f_sife, mtd.N, mat_b, mat_f);

    # compute the mean and mean diff of mat_b and mat_f
    bio.auxil.mat_mean .= (bio.auxil.mat_b .+ bio.auxil.mat_f) ./ 2;
    bio.auxil.mat_diff .= (bio.auxil.mat_b .- bio.auxil.mat_f) ./ 2;

    return nothing
);
