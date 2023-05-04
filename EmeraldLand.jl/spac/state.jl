#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to synchronize state among spac and spac state structs
#
#######################################################################################################################################################################################################
"""

    spac_state!(spac::MultiLayerSPAC{FT}, state::MultiLayerSPACState{FT}) where {FT<:AbstractFloat}
    spac_state!(state::MultiLayerSPACState{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Synchronize state variables from 1st to 2nd struct, given
- `spac` `MultiLayerSPAC` struct for SPAC
- `state` `MultiLayerSPACState` struct for states

"""
function spac_state! end

spac_state!(spac::MultiLayerSPAC{FT}, state::MultiLayerSPACState{FT}) where {FT<:AbstractFloat} = (
    (; CANOPY, DIM_LAYER, LEAVES, MEMORY) = spac;

    for _i in 1:DIM_LAYER
        state.gs_shaded[_i] = LEAVES[_i].g_H₂O_s_shaded;
        state.gs_sunlit[:,:,_i] .= LEAVES[_i].g_H₂O_s_sunlit;
    end;

    state.t_clm .= MEMORY.tem;

    # save the variables used for publications
    state.gpp = GPP(spac);
    state.tropomi_sif₆₈₃ = TROPOMI_SIF683(CANOPY);
    state.tropomi_sif₇₄₀ = TROPOMI_SIF740(CANOPY);

    return nothing
);

spac_state!(state::MultiLayerSPACState{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat} = (
    (; DIM_LAYER, LEAVES, MEMORY) = spac;

    for _i in 1:DIM_LAYER
        LEAVES[_i].g_H₂O_s_shaded = state.gs_shaded[_i];
        LEAVES[_i].g_H₂O_s_sunlit .= state.gs_sunlit[:,:,_i];
    end;

    MEMORY.tem .= state.t_clm;

    return nothing
);
