#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: rename the function from leaf_ETR! to to photosystem_electron_transport! to be more specific
#     2022-Jan-18: use ppar and p_i as inputs rather than field from Leaf, and this allows for more modular operations
#     2022-Jan-18: use ppar rather than par as in Johnson and Berry's paper (they convert par to ppar later)
#     2022-Feb-07: add e2c calculation
#     2022-Mar-01: save PSI J
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#
#######################################################################################################################################################################################################
"""

    photosystem_electron_transport!(psm::LeafPhotosystem{FT}, ppar::FT, p_i::FT; β::FT = FT(1)) where {FT}

Update the electron transport rates, given
- `psm` `LeafPhotosystem` type struct
- `ppar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`, used to compute e_to_c
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function photosystem_electron_transport! end;

photosystem_electron_transport!(psm::LeafPhotosystem{FT}, ppar::FT, p_i::FT; β::FT = FT(1)) where {FT} = photosystem_electron_transport!(psm.trait, psm.state, psm.auxil, ppar, p_i; β = β);

photosystem_electron_transport!(pst::C3CytoTrait{FT}, pss::C3CytoState{FT}, psa::PSMAuxil{FT}, ppar::FT, p_i::FT; β::FT = FT(1)) where {FT} = (
    psa.e2c   = (p_i - psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star);
    psa.j_psi = colimited_rate(β * psa.v_qmax, ppar * (1 - psa.f_psii) * psa.ϕ_psi_max, pst.COLIMIT_J);
    psa.η     = 1 - psa.η_l / psa.η_c + (3 * p_i + 7 * psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star) / psa.η_c;
    psa.j_pot = psa.j_psi / psa.η;

    return nothing
);

photosystem_electron_transport!(pst::C3VJPTrait{FT}, pss::C3VJPState{FT}, psa::PSMAuxil{FT}, ppar::FT, p_i::FT; β::FT = FT(1)) where {FT} = (
    psa.e2c   = (p_i - psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star);
    psa.j_pot = psa.f_psii * psa.ϕ_psii_max * ppar;
    psa.j     = colimited_rate(psa.j_pot, β * psa.j_max, pst.COLIMIT_J);

    return nothing
);

photosystem_electron_transport!(pst::Union{C4CLMTrait{FT}, C4VJPTrait{FT}}, pss::C4VJPState{FT}, psa::PSMAuxil{FT}, ppar::FT, p_i::FT; β::FT = FT(1)) where {FT} = (
    psa.e2c   = 1 / 6;
    psa.j_pot = psa.f_psii * psa.ϕ_psii_max * ppar;

    return nothing
);
