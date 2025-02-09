#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2025-Feb-09: add function to derive the SIF spectra for one plate using the doubling adding method
#
#######################################################################################################################################################################################################
"""

kubelka_munk_sif_matrices!(
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
function kubelka_munk_sif_matrices!(
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
        if ρ_no[i] + τ_no[i] <= 1
            d = sqrt((1 + ρ_no[i] + τ_no[i]) * (1 + ρ_no[i] - τ_no[i]) * (1 - ρ_no[i] + τ_no[i]) * (1 - ρ_no[i] - τ_no[i]));
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

    # 3. indices of WLE and WLF within wlp
    ϵ = FT(2) ^ -N;
    ρ_e     = s[IΛ_SIFE] .* ϵ;
    ρ_f     = s[IΛ_SIF] .* ϵ;
    sigmoid = 1 ./ (1 .+ exp.(-Λ_SIF ./ 10) .* exp.(Λ_SIFE' ./ 10));
    mat_f_  = Φ_PS[IΛ_SIF] .* ϵ ./ 2 .* k_chl[IΛ_SIFE]' .* sigmoid;
    mat_b_  = Φ_PS[IΛ_SIF] .* ϵ ./ 2 .* k_chl[IΛ_SIFE]' .* sigmoid;
    τ_e     = 1 .- (k[IΛ_SIFE] .+ s[IΛ_SIFE]) .* ϵ;
    τ_f     = 1 .- (k[IΛ_SIF] .+ s[IΛ_SIF]) .* ϵ;

    # 4. Doubling adding routine
    m_1_e = ones(FT, 1, length(IΛ_SIFE));
    m_f_1 = ones(FT, length(IΛ_SIF), 1);
    for _ in 1:N
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
        mat_f_n = mat_f_ .* a₁₁ .+ mat_b_ .* a₁₂;
        mat_b_n = mat_b_ .* a₂₁ .+ mat_f_ .* a₂₂;
        τ_e     = τ_e_n;
        ρ_e     = ρ_e_n;
        τ_f     = τ_f_n;
        ρ_f     = ρ_f_n;
        mat_f_  = mat_f_n;
        mat_b_  = mat_b_n;
    end;

    # 5. This reduced red SIF quite a bit in backscatter, not sure why.
    ρ_b  = ρ_no .+ τ_no .^ 2 .* ρ_i_21 ./ (1 .- ρ_no .* ρ_i_21);
    z_x  = τ_i_θ[IΛ_SIFE] ./ (1 .- ρ_i_21[IΛ_SIFE] .* ρ_b[IΛ_SIFE]);
    m_xe = m_f_1 .* z_x';
    m_xf = τ_i_21[IΛ_SIF] ./ (1 .- ρ_i_21[IΛ_SIF] .* ρ_b[IΛ_SIF]) .* m_1_e;
    z_y  = τ_no[IΛ_SIFE] .* ρ_i_21[IΛ_SIFE] ./ (1 .- ρ_no[IΛ_SIFE] .* ρ_i_21[IΛ_SIFE]);
    m_ye = m_f_1 .* z_y';
    m_yf = τ_no[IΛ_SIF] .* ρ_i_21[IΛ_SIF] ./ (1 .- ρ_no[IΛ_SIF] .* ρ_i_21[IΛ_SIF]) .* m_1_e;
    ma   = m_xe .* (1 .+ m_ye .* m_yf) .* m_xf;
    mb   = m_xe .* (m_ye .+ m_yf) .* m_xf;

    mat_b .= ma .* mat_b_ .+ mb .* mat_f_;
    mat_f .= ma .* mat_f_ .+ mb .* mat_b_;

    return nothing
end;
