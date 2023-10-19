module Optics

using ..Constant: AVOGADRO, H_PLANCK, LIGHT_SPEED


const FAC = 1e-9 / (H_PLANCK() * LIGHT_SPEED() * AVOGADRO());


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-22: rename the function to photon
#     2021-Oct-22: add a method to convert direct from number to number
#
#######################################################################################################################################################################################################
"""

    photon(λ::FT, E::FT) where {FT}

Return the number of moles of photons, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `E` Joules of energy

"""
function photon(λ::FT, E::FT) where {FT}
    return E * λ * FT(FAC)
end;


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

    photon!(λ::Vector{FT}, E::Vector{FT}, phot::Vector{FT}) where {FT}
    photon!(λ::Vector{FT}, E::Vector{FT}) where {FT}

Compute and save the number of moles of photons, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `E` Joules of energy (will be converted to moles of photons if phot in not given)
- `phot` Mole of photons (variable to save)

"""
function photon! end;

photon!(λ::Vector{FT}, E::Vector{FT}, phot::Vector{FT}) where {FT} = (phot .= photon.(λ, E); return nothing);

photon!(λ::Vector{FT}, E::Vector{FT}) where {FT} = (E .*= λ .* FT(FAC); return nothing);


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-22: define function to convert photon back to energy
#
#######################################################################################################################################################################################################
"""

    energy(λ::FT, phot::FT) where {FT}

Return the energy, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `phot` Number of moles of photon

"""
function energy(λ::FT, phot::FT) where {FT}
    return phot / (λ * FT(FAC))
end;


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

    energy!(λ::Vector{FT}, phot::Vector{FT}, E::Vector{FT}) where {FT}
    energy!(λ::Vector{FT}, phot::Vector{FT}) where {FT}

Compute and save the number of moles of photons, given
- `λ` Wave length in `[nm]`, converted to `[m]` by FAC
- `phot` Mole of photons (will be converted to moles of photons if E is not given)
- `E` Joules of energy (variable to save)

"""
function energy! end;

energy!(λ::Vector{FT}, phot::Vector{FT}, E::Vector{FT}) where {FT} = (E .= energy.(λ, phot); return nothing);

energy!(λ::Vector{FT}, phot::Vector{FT}) where {FT} = (phot ./= λ .* FT(FAC); return nothing);


end; # module
