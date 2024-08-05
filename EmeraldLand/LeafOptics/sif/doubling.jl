#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-21: add function to compute the matrices using doubling adding method
#
#######################################################################################################################################################################################################
"""

    doubling_sif_matrices(config::SPACConfiguration{FT}, bio::LeafBio{FT}; NDUB = 10) where {FT}

Return the matrices for fluorescence calculation using doubling adding method, given
- `config` `SPACConfiguration` type struct
- `bio` `LeafBio` type struct
- `NDUB` Number of doubling adding steps

"""
function doubling_sif_matrices(config::SPACConfiguration{FT}, bio::LeafBio{FT}; NDUB = 10) where {FT}
    (; SPECTRA) = config;
    (; IΛ_SIF, IΛ_SIFE, Λ_SIF, Λ_SIFE, Φ_PS) = SPECTRA;
    auxil = bio.auxil;

    # Doubling method used to calculate fluoresence is now only applied to the part of the leaf where absorption takes place, that is, the part exclusive of the leaf-air interfaces.
    # The reflectance (rho) and transmittance (tau) of this part of the leaf are now determined by "subtracting" the interfaces.
    # CF Note: All of the below takes about 10 times more time than the RT above. Need to rething speed and accuracy. (10nm is bringing it down a lot!)
    ρ_b = (auxil.ρ_leaf .- auxil.ρ_interface_θ) ./ (auxil.τ_interface_θ .* auxil.τ_interface_21 .+ (auxil.ρ_leaf - auxil.ρ_interface_θ) .* auxil.ρ_interface_21);
    z   = auxil.τ_leaf .* (1 .- ρ_b .* auxil.ρ_interface_21) ./ (auxil.τ_interface_θ .* auxil.τ_interface_21);
    ρ   = max.(0, (ρ_b .- auxil.ρ_interface_21 .* z .^ 2) ./ (1 .- (auxil.ρ_interface_21.* z) .^ 2));
    τ   = (1 .- ρ_b .* auxil.ρ_interface_21) ./ (1 .- (auxil.ρ_interface_21.* z) .^ 2) .* z;

    # Derive Kubelka-Munk s and k
    a = similar(ρ);
    b = similar(ρ);
    for i in eachindex(ρ)
        if ρ[i] + τ[i] < 1
            d = sqrt((1 + ρ[i] + τ[i]) * (1 + ρ[i] - τ[i]) * (1 - ρ[i] + τ[i]) *  (1 - ρ[i] - τ[i]));
            a[i] = (1 + ρ[i] ^ 2 - τ[i] ^ 2 + d) / (2 * ρ[i]);
            b[i] = (1 - ρ[i] ^ 2 + τ[i] ^ 2 + d) / (2 * τ[i]);
        else
            a[i] = 1;
            b[i] = 1;
        end;
    end;

    s = ρ ./ τ;
    k = log.(b);
    for i in eachindex(a)
        if 1 < a[i] < Inf
            s[i] = 2 * a[i] / (a[i] ^ 2 - 1) * log(b[i]);
            k[i] = (a[i] - 1) / (a[i] + 1) * log(b[i]);
        end;
    end;
    k_chl = auxil.f_sife .* k;

    # indices of WLE and WLF within wlp
    ϵ = FT(2) ^ -NDUB;
    ρ_e     = s[IΛ_SIFE] .* ϵ;
    ρ_f     = s[IΛ_SIF] .* ϵ;
    sigmoid = 1 ./ (1 .+ exp.(-Λ_SIF ./ 10) .* exp.(Λ_SIFE' ./ 10));
    mat_f   = Φ_PS[IΛ_SIF] .* ϵ ./ 2 .* k_chl[IΛ_SIFE]' .* sigmoid;
    mat_b   = Φ_PS[IΛ_SIF] .* ϵ ./ 2 .* k_chl[IΛ_SIFE]' .* sigmoid;
    τ_e     = 1 .- (k[IΛ_SIFE] .+ s[IΛ_SIFE]) .* ϵ;
    τ_f     = 1 .- (k[IΛ_SIF] .+ s[IΛ_SIF]) .* ϵ;

    # Doubling adding routine
    m_1_e = ones(FT, 1, length(IΛ_SIFE));
    m_f_1 = ones(FT, length(IΛ_SIF), 1);
    for _ in 1:NDUB
        x_e     = τ_e ./ (1 .- ρ_e .^ 2);
        x_f     = τ_f ./ (1 .- ρ_f .^ 2);
        τ_e_n   = τ_e .* x_e;
        τ_f_n   = τ_f .* x_f;
        ρ_e_n   = ρ_e .* (1 .+ τ_e_n);
        ρ_f_n   = ρ_f .* (1 .+ τ_f_n);
        a₁₁     = x_f .* m_1_e .+ m_f_1 .* x_e';
        a₁₂     = (x_f .* x_e') .* (ρ_f .* m_1_e .+ m_f_1 .* ρ_e');
        a₂₁     = 1 .+ (x_f * x_e') .* (1 .+ ρ_f * ρ_e');
        z_e     = x_e .* ρ_e;
        z_f     = x_f .* ρ_f;
        a₂₂     = z_f .* m_1_e .+ m_f_1 .* z_e';
        mat_f_n = mat_f .* a₁₁ .+ mat_b .* a₁₂;
        mat_b_n = mat_b .* a₂₁ .+ mat_f .* a₂₂;
        τ_e     = τ_e_n;
        ρ_e     = ρ_e_n;
        τ_f     = τ_f_n;
        ρ_f     = ρ_f_n;
        mat_f   = mat_f_n;
        mat_b   = mat_b_n;
    end;

    # This reduced red SIF quite a bit in backscatter, not sure why.
    ρ_b  = ρ .+ τ .^ 2 .* auxil.ρ_interface_21 ./ (1 .- ρ .* auxil.ρ_interface_21);
    z_x  = auxil.τ_interface_θ[IΛ_SIFE] ./ (1 .- auxil.ρ_interface_21[IΛ_SIFE] .* ρ_b[IΛ_SIFE]);
    m_xe = m_f_1 .* z_x';
    m_xf = auxil.τ_interface_21[IΛ_SIF] ./ (1 .- auxil.ρ_interface_21[IΛ_SIF] .* ρ_b[IΛ_SIF]) .* m_1_e;
    z_y  = τ[IΛ_SIFE] .* auxil.ρ_interface_21[IΛ_SIFE] ./ (1 .- ρ[IΛ_SIFE] .* auxil.ρ_interface_21[IΛ_SIFE]);
    m_ye = m_f_1 .* z_y';
    m_yf = τ[IΛ_SIF] .* auxil.ρ_interface_21[IΛ_SIF] ./ (1 .- ρ[IΛ_SIF] .* auxil.ρ_interface_21[IΛ_SIF]) .* m_1_e;
    ma   = m_xe .* (1 .+ m_ye .* m_yf) .* m_xf;
    mb   = m_xe .* (m_ye .+ m_yf) .* m_xf;

    leaf_mat_b = ma .* mat_b .+ mb .* mat_f;
    leaf_mat_f = ma .* mat_f .+ mb .* mat_b;

    return leaf_mat_b, leaf_mat_f
end;
