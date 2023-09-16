#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute reflectance and transmittance at the interface for direct light
#
#######################################################################################################################################################################################################
"""

    interface_ρ_τ(n₁::FT, n₂::FT, θ₁::FT)

Return the reflectance and transmittance at the interface between two media for direct radiation, given
- `n₁` refractive index of the incoming radiation medium
- `n₂` refractive index of the outgoing radiation medium
- `θ₁` incident angle of the incoming radiation

"""
function interface_ρ_τ(n₁::FT, n₂::FT, θ₁::FT) where {FT}
    sinθ₁ = sind(θ₁);
    cosθ₁ = cosd(θ₁);
    sinθ₂ = min(1, n₁ / n₂ * sinθ₁);
    cosθ₂ = sqrt(1 - sinθ₂^2);

    ρ_s = ((n₁ * cosθ₁ - n₂ * cosθ₂) / (n₁ * cosθ₁ + n₂ * cosθ₂)) ^ 2;
    ρ_p = ((n₂ * cosθ₁ - n₁ * cosθ₂) / (n₂ * cosθ₁ + n₁ * cosθ₂)) ^ 2;
    ρ = (ρ_s + ρ_p) / 2;

    return ρ, 1 - ρ
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute transmittance at the interface for light integrated from 0 to a certain degree
#
#######################################################################################################################################################################################################
"""

    interface_integrated_ρ_τ(n₁::FT, n₂::FT, θ₁::FT)

Return the integrated transmittance (two polarations) at the interface between two media for light integrated from 0 to a certain degree, given
- `m` refractive index ratio of the two media
- `θ₁` incident angle of the incoming radiation

"""
function interface_integrated_τ(m::FT, θ₁::FT) where {FT}
    @assert m > 1;

    m² = m ^ 2;
    m⁴ = m ^ 4;
    sinθ₁ = sind(θ₁);
    cosθ₁ = cosd(θ₁);
    sinθ₂ = sinθ₁ / m;
    cosθ₂ = sqrt(1 - sinθ₂^2);

    sin²θ₁ = sinθ₁ ^ 2;
    sin⁴θ₁ = sinθ₁ ^ 4;
    u = 1 + 2 * cosθ₁ * (cosθ₁ + m * cosθ₂) / (m² - 1);
    U = (m² + 1) * u - (m² - 1);

    s₀ = FT(-2 / 3) * FT(π) * (m² - 1) ^ -2;
    s₁ = (2 * m * cosθ₁ * cosθ₂) ^ 3;
    s₂ = 4sin²θ₁ * (2sin⁴θ₁ - (3m² + 3) * sin²θ₁ + 6m²);
    i_s = s₀ * (s₁ + s₂);

    p₁ = 2 * FT(π) * m² * (m² - 1) * (m² + 1) ^ -3;
    p₂ = U + (m⁴ + 6m² + 1) / U + 2 * (m² - 1) * log(U);
    p₃ = 2 * FT(π) * m² * (m² - 1) ^ -2;
    p₄ = U / u - (m² + 1)^2 * u / U - 2 * (m² + 1) * log(U/u);
    i_p = p₁ * p₂ + p₃ * p₄;

    return i_s, i_p
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute transmittance at the interface for isotropic light
#     2023-Sep-15: fix a typo introduced by mistake when removing _ from the variable names
#
#######################################################################################################################################################################################################
"""

    interface_isotropic_τ(n₁::FT, n₂::FT, θ₁::FT)

Return the transmittance at the interface between two media for isotropic radiation, given
- `n₁` refractive index of the incoming radiation medium
- `n₂` refractive index of the outgoing radiation medium
- `θ₁` incident angle of the incoming radiation

"""
function interface_isotropic_τ(n₁::FT, n₂::FT, θ₁::FT) where {FT}
    m = n₂ / n₁;

    if m > 1
        s₁, p₁ = interface_integrated_τ(m, FT(0));
        s₂, p₂ = interface_integrated_τ(m, θ₁);
        _denom = FT(2) * FT(π) * sind(θ₁)^2;

        τ_s = (s₂ - s₁) / _denom;
        τ_p = (p₁ - p₂) / _denom; # seems like the original equation made a mistake in the sign

        return (τ_s + τ_p) / 2
    else
        s₁, p₁ = interface_integrated_τ(1 / m, FT(0));
        s₂, p₂ = interface_integrated_τ(1 / m, θ₁);
        _denom = FT(2) * FT(π) * sind(θ₁)^2;

        τ_s = (s₂ - s₁) / _denom;
        τ_p = (p₁ - p₂) / _denom;

        return (τ_s + τ_p) / 2 * m ^ 2
    end;
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to update the auxiliary variables in HyperLeafBio struct
#
#######################################################################################################################################################################################################
"""

    leaf_interface_ρ_τ!(config::SPACConfiguration{FT}, bio::HyperLeafBio{FT}, θ::FT) where {FT}
    leaf_interface_ρ_τ!(lha::HyperspectralAbsorption{FT}, bio::HyperLeafBio{FT}, θ::FT) where {FT}

Update the interface reflectance and transmittance in `bio`, given
- `config` SPAC configuration
- `bio` HyperLeafBio struct
- `θ` average angle of the incident radiation
- `lha` HyperspectralAbsorption struct

"""
function leaf_interface_ρ_τ! end;

leaf_interface_ρ_τ!(config::SPACConfiguration{FT}, bio::HyperLeafBio{FT}, θ::FT) where {FT} = leaf_interface_ρ_τ!(config.LHA, bio, θ);

leaf_interface_ρ_τ!(lha::HyperspectralAbsorption{FT}, bio::HyperLeafBio{FT}, θ::FT) where {FT} = (
    (; NR) = lha;

    bio.auxil.τ_interface_θ  .= interface_isotropic_τ.(FT(1), NR, θ);
    bio.auxil.τ_interface_12 .= interface_isotropic_τ.(FT(1), NR, FT(90));
    bio.auxil.τ_interface_21 .= interface_isotropic_τ.(NR, FT(1), FT(90));

    bio.auxil.ρ_interface_θ  .= 1 .- bio.auxil.τ_interface_θ;
    bio.auxil.ρ_interface_12 .= 1 .- bio.auxil.τ_interface_12;
    bio.auxil.ρ_interface_21 .= 1 .- bio.auxil.τ_interface_21;

    return nothing
);
