module Optics

using ..Constant: AVOGADRO, H_PLANCK, LIGHT_SPEED


const FAC = 1e-9 / (H_PLANCK() * LIGHT_SPEED() * AVOGADRO());


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2020-Mar-30: migrate the function from SCOPE
#     2021-Oct-21: rename the function to average_transmittance
# To do
#     TODO: figure out where does this equation comes from
#
#######################################################################################################################################################################################################
"""

    average_transmittance(α::FT, nr::FT) where {FT<:AbstractFloat}

Return the average transmittance of isotropic radiation across an interface between two dielectrics, given
- `α` angle of incidence
- `nr` Index of refraction

# References
- Stern (1964) Transmission of isotropic radiation across an interface between two dielectrics. Applied Optics 3(1): 111-113.
- Allen (1973) Transmission of isotropic light across a dielectric surface in two and three dimensions. Journal of the Optical Society of America 63(6): 664-666.

"""
function average_transmittance(α::FT, nr::FT) where {FT<:AbstractFloat}
    @assert 0 < α <= 90;

    # some shortcuts to avoid overly comlicated equation
    _a     = (nr + 1) ^ 2 / 2;
    _a³    = _a ^ 3;
    _n²    = nr ^ 2;
    _n⁴    = nr ^ 4;
    _n⁶    = nr ^ 6;
    _n²p   = _n² + 1;
    _n²p²  = _n²p ^ 2;
    _n²p³  = _n²p ^ 3;
    _n²m²  = (_n² - 1) ^ 2;
    _k     = -1 * _n²m² / 4;
    _k²    = _k ^ 2;
    _sin²α = sind(α) ^ 2;
    _b₂    = _sin²α - _n²p/2;
    _b₁    = (α==90 ? 0 : sqrt(_b₂ ^ 2 + _k));
    _b     = _b₁ - _b₂;
    _b³    = _b ^ 3;
    _npanm = 2 * _n²p * _a - _n²m²;
    _npbnm = 2 * _n²p * _b - _n²m²;

    # S polarization
    _ts  = ( _k² / (6*_b³) + _k/_b - _b/2 ) - ( _k² / (6*_a³) + _k/_a - _a/2 );

    # P polarization
    _tp₁ = -2 * _n² * (_b - _a) / _n²p²;
    _tp₂ = -2 * _n² * _n²p * log(_b / _a) / _n²m²;
    _tp₃ = _n² * (1/_b - 1/_a) / 2;
    _tp₄ = 16 * _n⁴ * (_n⁴ + 1) * log(_npbnm / _npanm) / (_n²p³ * _n²m²);
    _tp₅ = 16 * _n⁶ * (1/_npbnm - 1/_npanm) / _n²p³;
    _tp  = _tp₁ + _tp₂ + _tp₃ + _tp₄ + _tp₅;

    return  (_ts + _tp) / (2 * _sin²α)
end


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-22: rename the function to photon
#     2021-Oct-22: add a method to convert direct from number to number
#
#######################################################################################################################################################################################################
"""

    photon(λ::FT, E::FT) where {FT<:AbstractFloat}

Return the number of moles of photons, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `E` Joules of energy

"""
function photon(λ::FT, E::FT) where {FT<:AbstractFloat}
    return E * λ * FT(FAC)
end


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2022-Jun-13: add function
#     2021-Jun-13: add method to save to provided 3rd variable
#     2021-Jun-13: add method to save to provided 2rd variable
#
#######################################################################################################################################################################################################
"""

    photon!(λ::Vector{FT}, E::Vector{FT}, phot::Vector{FT}) where {FT<:AbstractFloat}
    photon!(λ::Vector{FT}, E::Vector{FT}) where {FT<:AbstractFloat}

Compute and save the number of moles of photons, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `E` Joules of energy (will be converted to moles of photons if phot in not given)
- `phot` Mole of photons (variable to save)

"""
function photon! end

photon!(λ::Vector{FT}, E::Vector{FT}, phot::Vector{FT}) where {FT<:AbstractFloat} = (phot .= photon.(λ, E); return nothing);

photon!(λ::Vector{FT}, E::Vector{FT}) where {FT<:AbstractFloat} = (E .*= λ .* FT(FAC); return nothing);


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-22: define function to convert photon back to energy
#
#######################################################################################################################################################################################################
"""

    energy(λ::FT, phot::FT) where {FT<:AbstractFloat}

Return the energy, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `phot` Number of moles of photon

"""
function energy(λ::FT, phot::FT) where {FT<:AbstractFloat}
    return phot / (λ * FT(FAC))
end


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2022-Jun-13: add function
#     2021-Jun-13: add method to save to provided 3rd variable
#     2021-Jun-13: add method to save to provided 2rd variable
#
#######################################################################################################################################################################################################
"""

    energy!(λ::Vector{FT}, phot::Vector{FT}, E::Vector{FT}) where {FT<:AbstractFloat}
    energy!(λ::Vector{FT}, phot::Vector{FT}) where {FT<:AbstractFloat}

Compute and save the number of moles of photons, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `phot` Mole of photons (will be converted to moles of photons if E is not given)
- `E` Joules of energy (variable to save)

"""
function energy! end

energy!(λ::Vector{FT}, phot::Vector{FT}, E::Vector{FT}) where {FT<:AbstractFloat} = (E .= energy.(λ, phot); return nothing);

energy!(λ::Vector{FT}, phot::Vector{FT}) where {FT<:AbstractFloat} = (phot ./= λ .* FT(FAC); return nothing);


end # module
