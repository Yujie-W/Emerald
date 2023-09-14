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
    _sin_θ₁ = sind(θ₁);
    _cos_θ₁ = cosd(θ₁);
    _sin_θ₂ = min(1, n₁ / n₂ * _sin_θ₁);
    _cos_θ₂ = sqrt(1 - _sin_θ₂^2);

    _ρ_s = ((n₁ * _cos_θ₁ - n₂ * _cos_θ₂) / (n₁ * _cos_θ₁ + n₂ * _cos_θ₂)) ^ 2;
    _ρ_p = ((n₂ * _cos_θ₁ - n₁ * _cos_θ₂) / (n₂ * _cos_θ₁ + n₁ * _cos_θ₂)) ^ 2;
    _ρ = (_ρ_s + _ρ_p) / 2;

    return _ρ, 1 - _ρ
end


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

    _s_0 = FT(-2 / 3) * FT(π) * (m² - 1) ^ -2;
    _s_1 = (2 * m * cosθ₁ * cosθ₂) ^ 3;
    _s_2 = 4sin²θ₁ * (2sin⁴θ₁ - (3m² + 3) * sin²θ₁ + 6m²);
    _i_s = _s_0 * (_s_1 + _s_2);

    _p_1 = 2 * FT(π) * m² * (m² - 1) * (m² + 1) ^ -3;
    _p_2 = U + (m⁴ + 6m² + 1) / U + 2 * (m² - 1) * log(U);
    _p_3 = 2 * FT(π) * m² * (m² - 1) ^ -2;
    _p_4 = U / u - (m² + 1)^2 * u / U - 2 * (m² + 1) * log(U/u);
    _i_p = _p_1 * _p_2 + _p_3 * _p_4;

    return _i_s, _i_p
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-14: add function to compute transmittance at the interface for isotropic light
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
        _s_1, _p_1 = interface_integrated_τ(m, FT(0));
        _s_2, _p_2 = interface_integrated_τ(m, θ₁);
        _denom = FT(2) * FT(π) * sind(θ₁)^2;

        _τ_s = (_s_2 - _s_1) / _denom;
        _τ_p = (_p_1 - _p_2) / _denom; # seems like the original equation made a mistake in the sign
        _τ = (_τ_s + _τ_p) / 2;
    else
        _s_1, _p_1 = interface_integrated_τ(1 / m, FT(0));
        _s_2, _p_2 = interface_integrated_τ(1 / m, θ₁);
        _denom = FT(2) * FT(π) * sind(θ₁)^2;

        _τ_s = (_s_2 - _s_1) / _denom;
        _τ_p = (_p_1 - _p_2) / _denom;
        _τ = (_τ_s + _τ_p) / 2 * m ^ 2;
    end

    return _τ
end
