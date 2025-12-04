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
#     2024-Jul-22: save j as j_pot fpr C4, C3Cyto models
#     2024-Aug-01: generalize the function for GeneralC3Trait and GeneralC4Trait
#     2024-Aug-05: make sure COLIMIT_J is SerialColimit for AjMethodC3VqmaxPi
#
#######################################################################################################################################################################################################
"""

    photosystem_electron_transport!(
                cache::SPACCache{FT},
                ps::CanopyLayerPhotosystem{FT},
                ppar::Vector{FT},
                p_i::Union{FT, Vector{FT}};
                β::FT = FT(1)) where {FT}
    photosystem_electron_transport!(
                ps::LeafPhotosystem{FT},
                ppar::FT,
                p_i::FT;
                β::FT = FT(1)) where {FT}

Update the electron transport rates, given
- `cache` `SPACCache` type struct
- `ps` `CanopyLayerPhotosystem` or `LeafPhotosystem` type struct
- `ppar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`, used to compute e_to_c
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function photosystem_electron_transport! end;

# For Leaf
photosystem_electron_transport!(
            ps::LeafPhotosystem{FT},
            ppar::FT,
            p_i::FT;
            β::FT = FT(1)) where {FT} = photosystem_electron_transport!(ps.trait, ps.state, ps.auxil, ppar, p_i; β = β);

photosystem_electron_transport!(
            pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::LeafPhotosystemAuxil{FT},
            ppar::FT,
            p_i::FT;
            β::FT = FT(1)) where {FT} = photosystem_electron_transport!(pst, pss, psa, pst.AJM, ppar, p_i; β = β);

photosystem_electron_transport!(
            pst::GeneralC3Trait{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ajm::AjMethodC3JmaxPi,
            ppar::FT,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (
    psa.e2c   = (p_i - psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star);
    psa.j_pot = psa.f_psii * psa.ϕ_psii_max * ppar;
    psa.j     = colimited_rate(psa.j_pot, β * psa.j_max, pst.COLIMIT_J);

    return nothing
);

photosystem_electron_transport!(
            pst::GeneralC3Trait{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ajm::AjMethodC3VqmaxPi,
            ppar::FT,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (
    @assert pst.COLIMIT_J isa SerialColimit "J Limitation must be serial colimit";

    psa.e2c   = (p_i - psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star);
    psa.j_psi = colimited_rate(β * psa.v_qmax, ppar * (1 - psa.f_psii) * psa.ϕ_psi_max, pst.COLIMIT_J);
    psa.η     = 1 - psa.η_l / psa.η_c + (3 * p_i + 7 * psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star) / psa.η_c;
    psa.j_pot = psa.j_psi / psa.η;
    psa.j     = psa.j_pot;

    return nothing
);

photosystem_electron_transport!(
            pst::GeneralC4Trait{FT},
            pss::C4State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ajm::AjMethodC4JPSII,
            ppar::FT,
            p_i::FT;
            β::FT = FT(1)) where {FT} = (
    psa.e2c   = 1 / 6;
    psa.j_pot = psa.f_psii * psa.ϕ_psii_max * ppar;
    psa.j     = psa.j_pot;

    return nothing
);

# For CanopyLayer
photosystem_electron_transport!(
            cache::SPACCache{FT},
            ps::CanopyLayerPhotosystem{FT},
            ppar::Vector{FT},
            p_i::Union{FT, Vector{FT}};
            β::FT = FT(1)) where {FT} = photosystem_electron_transport!(cache, ps.trait, ps.state, ps.auxil, ppar, p_i; β = β);

photosystem_electron_transport!(
            cache::SPACCache{FT},
            pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ppar::Vector{FT},
            p_i::Union{FT, Vector{FT}};
            β::FT = FT(1)) where {FT} = photosystem_electron_transport!(cache, pst, pss, psa, pst.AJM, ppar, p_i; β = β);

photosystem_electron_transport!(
            cache::SPACCache{FT},
            pst::GeneralC3Trait{FT},
            pss::C3State{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ajm::AjMethodC3JmaxPi,
            ppar::Vector{FT},
            p_i::Union{FT, Vector{FT}};
            β::FT = FT(1)) where {FT} = (
    @. psa.e2c   = (p_i - psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star);
    @. psa.j_pot = psa.f_psii * psa.ϕ_psii_max * ppar;
    colimited_rate!(β * psa.j_max, psa.j_pot, psa.j, pst.COLIMIT_J);

    return nothing
);

photosystem_electron_transport!(
            cache::SPACCache{FT},
            pst::GeneralC3Trait{FT},
            pss::C3State{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ajm::AjMethodC3VqmaxPi,
            ppar::Vector{FT},
            p_i::Union{FT, Vector{FT}};
            β::FT = FT(1)) where {FT} = (
    @assert pst.COLIMIT_J isa SerialColimit "J Limitation must be serial colimit";

    _j = cache.cache_incl_azi_2_1;
    @. _j = ppar * (1 - psa.f_psii) * psa.ϕ_psi_max;
    colimited_rate!(β * psa.v_qmax, _j, psa.j_psi, pst.COLIMIT_J);

    @. psa.η     = 1 - psa.η_l / psa.η_c + (3 * p_i + 7 * psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star) / psa.η_c;
    @. psa.j_pot = psa.j_psi / psa.η;
    @. psa.j     = psa.j_pot;
    @. psa.e2c   = (p_i - psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star);

    return nothing
);

photosystem_electron_transport!(
            cache::SPACCache{FT},
            pst::GeneralC4Trait{FT},
            pss::C4State{FT},
            psa::CanopyLayerPhotosystemAuxil{FT},
            ajm::AjMethodC4JPSII,
            ppar::Vector{FT},
            p_i::Union{FT, Vector{FT}};
            β::FT = FT(1)) where {FT} = (
    @. psa.e2c   = 1 / 6;
    @. psa.j_pot = psa.f_psii * psa.ϕ_psii_max * ppar;
    @. psa.j     = psa.j_pot;

    return nothing
);
