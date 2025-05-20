#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2025-May-20: refactor the function with my own understanding, and found a bug in the original code...
#
#######################################################################################################################################################################################################
"""

    kubelka_munk_sif_matrices_new!(
                config::SPACConfiguration{FT},
                ρ_plate::Vector{FT},
                τ_plate::Vector{FT},
                ρ_i_θ::Vector{FT},
                τ_i_θ::Vector{FT},
                ρ_i_21::Vector{FT},
                τ_i_21::Vector{FT},
                f_sife::Vector{FT},
                N::Int,
                mat_b::Matrix{FT},
                mat_f::Matrix{FT}) where {FT}

Update the SIF matrices for one plate (or effective plate) using the doubling adding method, given
- `config` SPAC configuration
- `ρ_plate` reflectance of the plate
- `τ_plate` transmittance of the plate
- `ρ_i_θ` reflectance of the interface at incoming radiation angle (theta or diffuse)
- `τ_1_θ` transmittance of the interface at incoming radiation angle (theta or diffuse)
- `ρ_i_21` reflectance of the regular water-air interface
- `τ_i_21` transmittance of the regular water-air interface
- `f_sife` SIF excitation efficiency
- `N` number of doubling adding steps
- `mat_b` backward SIF matrix to be updated
- `mat_f` forward SIF matrix to be updated

"""
function kubelka_munk_sif_matrices_new!(
            config::SPACConfiguration{FT},
            ρ_plate::Vector{FT},
            τ_plate::Vector{FT},
            ρ_i_θ::Vector{FT},
            τ_i_θ::Vector{FT},
            ρ_i_21::Vector{FT},
            τ_i_21::Vector{FT},
            f_sife::Vector{FT},
            N::Int,
            mat_b::Matrix{FT},
            mat_f::Matrix{FT}) where {FT}
    (; IΛ_SIF, IΛ_SIFE, Λ_SIF, Λ_SIFE, Φ_PS) = config.SPECTRA;

    # 1. Get the rho and tau for the case of removing the top air-water interface (background reflectance)
    R_b  = @. (ρ_plate - ρ_i_θ) / (τ_i_θ * τ_i_21 + (ρ_plate - ρ_i_θ) * ρ_i_21);
    z    = @. τ_plate * (1 - R_b * ρ_i_21) / (τ_i_θ * τ_i_21);
    ρ_no = @. max(0, (R_b - ρ_i_21 * z ^ 2) / (1 - (ρ_i_21 * z) ^ 2));
    τ_no = @. (1 - R_b * ρ_i_21) / (1 - (ρ_i_21 * z) ^ 2) * z;

    # make sure the combined values are okay
    a = ones(FT, length(ρ_no));
    b = ones(FT, length(ρ_no));
    for i in eachindex(a)
        d² = (1 + ρ_no[i] + τ_no[i]) * (1 + ρ_no[i] - τ_no[i]) * (1 - ρ_no[i] + τ_no[i]) * (1 - ρ_no[i] - τ_no[i]);
        if d² > 0
            d = sqrt(d²);
            a[i] = (1 + ρ_no[i] ^ 2 - τ_no[i] ^ 2 + d) / (2 * ρ_no[i]);
            b[i] = (1 - ρ_no[i] ^ 2 + τ_no[i] ^ 2 + d) / (2 * τ_no[i]);
        end;
    end;

    # 2. Derive Kubelka-Munk s and k (as well as m from s and k)
    s = @. ρ_no / τ_no;
    k = @. log(b);
    for i in eachindex(s)
        if 1 < a[i] < Inf
            s[i] = 2 * a[i] / (a[i] ^ 2 - 1) * log(b[i]);
            k[i] = (a[i] - 1) / (a[i] + 1) * log(b[i]);
        end;
    end;
    k_chl = @. f_sife * k;

    # 3. the doubling rountine using for loop to loop through the excitation wavelength
    #    i for n-1  (start from 1)
    #    j for n (start from 2)
    f_b_i = ones(FT, length(IΛ_SIF));
    f_f_i = ones(FT, length(IΛ_SIF));
    f_b_j = ones(FT, length(IΛ_SIF));
    f_f_j = ones(FT, length(IΛ_SIF));
    ρ_f_i = ones(FT, length(IΛ_SIF));
    τ_f_i = ones(FT, length(IΛ_SIF));
    ρ_f_j = ones(FT, length(IΛ_SIF));
    τ_f_j = ones(FT, length(IΛ_SIF));
    ρ_e_i = FT(0);
    τ_e_i = FT(0);
    ρ_e_j = FT(0);
    τ_e_j = FT(0);

    ρ_i_21_f = ρ_i_21[IΛ_SIF];
    τ_i_21_f = τ_i_21[IΛ_SIF];
    τ_no_f = τ_no[IΛ_SIF];
    R_b_f = R_b[IΛ_SIF];

    δx = FT(2) ^ -N;
    ρ_f_1 = s[IΛ_SIF] .* δx;
    τ_f_1 = 1 .- (k[IΛ_SIF] .+ s[IΛ_SIF]) .* δx;
    for i in eachindex(IΛ_SIFE)
        i_e = IΛ_SIFE[i];
        ρ_e_1 = s[i_e] * δx;
        τ_e_1 = 1 - (k[i_e] + s[i_e]) * δx;
        sigmoid = @. 1 / (1 + exp(-Λ_SIF / 10) * exp(Λ_SIFE[i_e]/10));
        f_b_1 = Φ_PS[IΛ_SIF] ./ 2 .* k_chl[i_e] .* δx .* sigmoid;
        f_f_1 = Φ_PS[IΛ_SIF] ./ 2 .* k_chl[i_e] .* δx .* sigmoid;
        # initialization
        f_b_i .= f_b_1;
        f_f_i .= f_f_1;
        ρ_f_i .= ρ_f_1;
        τ_f_i .= τ_f_1;
        ρ_e_i  = ρ_e_1;
        τ_e_i  = τ_e_1;
        for _ in 1:N
            # for excitation doubling
            ρ_e_j = ρ_e_i + τ_e_i * ρ_e_i / (1 - ρ_e_i^2);
            τ_e_j = τ_e_i^2 / (1 - ρ_e_i^2);
            # for fluorescence doubling
            ρ_f_j = @. ρ_f_i + τ_f_i * ρ_f_i / (1 - ρ_f_i^2);
            τ_f_j = @. τ_f_i^2 / (1 - ρ_f_i^2);
            # for fluorescence scattering doubling
            f_b_j = @. f_b_i + f_f_i * ρ_f_i * τ_f_i / (1 - ρ_f_i^2) +
                       τ_e_i / (1 - ρ_e_i^2) * f_b_i * τ_f_i / (1 - ρ_f_i^2) +
                       τ_e_i * ρ_e_i / (1 - ρ_e_i^2) * (f_f_i + f_b_i * ρ_f_i * τ_f_i / (1 - ρ_f_i^2));
            f_f_j = @. f_f_i * τ_f_i / (1 - ρ_f_i^2) +
                       τ_e_i / (1 - ρ_e_i^2) * (f_f_i + f_b_i * ρ_f_i * τ_f_i / (1 - ρ_f_i^2)) +
                       τ_e_i * ρ_e_i / (1 - ρ_e_i^2) * f_b_i * τ_f_i / (1 - ρ_f_i^2);
            # update the matrices to the next level
            ρ_e_i = ρ_e_j;
            τ_e_i = τ_e_j;
            ρ_f_i = @. ρ_f_j;
            τ_f_i = @. τ_f_j;
            f_b_i = @. f_b_j;
            f_f_i = @. f_f_j;
        end;
        # account for the incoming radiation into the top of the plate without interface, and get the SIF without rescattering by the interfaces
        f_b_no = @. τ_i_θ[i_e] / (1 - ρ_i_21[i_e] * R_b[i_e]) * f_b_j;
        f_f_no = @. τ_i_θ[i_e] / (1 - ρ_i_21[i_e] * R_b[i_e]) * f_f_j;
        # account for the scattering by the interfaces
        f_b = @. (f_b_no + f_f_no * ρ_i_21_f * τ_no_f / (1 - ρ_i_21_f * τ_no_f)) * τ_i_21_f / (1 - ρ_i_21_f * R_b_f);
        f_f = @. (f_f_no + f_b_no * ρ_i_21_f * τ_no_f / (1 - ρ_i_21_f * τ_no_f)) * τ_i_21_f / (1 - ρ_i_21_f * R_b_f);
        # update the matrices
        mat_b[:,i] .= f_b;
        mat_f[:,i] .= f_f;
    end;

    return nothing
end;


test_leaf_sif_matrices!(config::SPACConfiguration{FT}, bio::LeafBio{FT}, mtd::SIFMatrixFluspectMethod) where {FT} = (
    (; ρ_leaf, τ_leaf, ρ_interface_θ, τ_interface_θ, ρ_interface_21, τ_interface_21, f_sife, mat_b, mat_f) = bio.auxil;

    # update the mat_b and mat_f based on the doubling adding method (which has been generalized to work for both Fluspect and Dualspect)
    kubelka_munk_sif_matrices!(config, ρ_leaf, τ_leaf, ρ_interface_θ, τ_interface_θ, ρ_interface_21, τ_interface_21, f_sife, mtd.N, mat_b, mat_f);
    mat_bb = deepcopy(mat_b);
    mat_ff = deepcopy(mat_f);
    kubelka_munk_sif_matrices_new!(config, ρ_leaf, τ_leaf, ρ_interface_θ, τ_interface_θ, ρ_interface_21, τ_interface_21, f_sife, mtd.N, mat_b, mat_f);

    return mat_bb, mat_ff, mat_b, mat_f
);
