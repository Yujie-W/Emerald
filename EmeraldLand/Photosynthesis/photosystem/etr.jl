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

    photosystem_electron_transport!(psm::C3Cyto{FT}, ppar::FT, p_i::FT; β::FT = FT(1)) where {FT}
    photosystem_electron_transport!(psm::C3VJP{FT}, ppar::FT, p_i::FT; β::FT = FT(1)) where {FT}
    photosystem_electron_transport!(psm::C4VJP{FT}, ppar::FT, p_i::FT; β::FT = FT(1)) where {FT}

Update the electron transport rates, given
- `psm` `C3Cyto`, `C3VJP`, or `C4VJP` type struct
- `ppar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`, used to compute e_to_c
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function photosystem_electron_transport! end;

photosystem_electron_transport!(psm::C3Cyto{FT}, ppar::FT, p_i::FT; β::FT = FT(1)) where {FT} = (
    psm.auxil.e2c   = (p_i - psm.auxil.γ_star) / (psm.state.EFF_1 * p_i + psm.state.EFF_2 * psm.auxil.γ_star);
    psm.auxil.j_psi = colimited_rate(β * psm.auxil.v_qmax, ppar * psm.state.F_PSI * psm.state.Φ_PSI_MAX, psm.state.COLIMIT_J);
    psm.auxil.η     = 1 - psm.auxil.η_l / psm.auxil.η_c + (3 * p_i + 7 * psm.auxil.γ_star) / (psm.state.EFF_1 * p_i + psm.state.EFF_2 * psm.auxil.γ_star) / psm.auxil.η_c;
    psm.auxil.j_pot = psm.auxil.j_psi / psm.auxil.η;

    return nothing
);

photosystem_electron_transport!(psm::C3VJP{FT}, ppar::FT, p_i::FT; β::FT = FT(1)) where {FT} = (
    psm.auxil.e2c   = (p_i - psm.auxil.γ_star) / (psm.state.EFF_1 * p_i + psm.state.EFF_2 * psm.auxil.γ_star);
    psm.auxil.j_pot = psm.state.F_PSII * psm.auxil.ϕ_psii_max * ppar;
    psm.auxil.j     = colimited_rate(psm.auxil.j_pot, β * psm.auxil.j_max, psm.state.COLIMIT_J);

    return nothing
);

photosystem_electron_transport!(psm::C4VJP{FT}, ppar::FT, p_i::FT; β::FT = FT(1)) where {FT} = (
    psm.auxil.e2c   = 1 / 6;
    psm.auxil.j_pot = psm.state.F_PSII * psm.auxil.ϕ_psii_max * ppar;

    return nothing
);
