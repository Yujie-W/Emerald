#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to synchronize state among spac and spac state structs
#     2023-Apr-13: add config to function call
#
#######################################################################################################################################################################################################
"""

    spac_state!(spac::MultiLayerSPAC{FT}, state::MultiLayerSPACState{FT}) where {FT}
    spac_state!(state::MultiLayerSPACState{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Synchronize state variables from 1st to 2nd struct, given
- `spac` `MultiLayerSPAC` struct for SPAC
- `state` `MultiLayerSPACState` struct for states

"""
function spac_state! end

spac_state!(spac::MultiLayerSPAC{FT}, state::MultiLayerSPACState{FT}, config::SPACConfiguration{FT}) where {FT} = (
    (; CANOPY, DIM_LAYER, LEAVES, MEMORY) = spac;
    (; WLSET) = config;

    for _i in 1:DIM_LAYER
        state.gs_shaded[_i] = LEAVES[_i].g_H₂O_s_shaded;
        state.gs_sunlit[:,:,_i] .= LEAVES[_i].g_H₂O_s_sunlit;
    end;

    state.t_clm .= MEMORY.tem;

    # save the variables used for publications
    state.gpp = GPP(spac);
    state.tropomi_sif₆₈₃ = TROPOMI_SIF683(CANOPY, WLSET);
    state.tropomi_sif₇₄₀ = TROPOMI_SIF740(CANOPY, WLSET);

    return nothing
);

spac_state!(state::MultiLayerSPACState{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; DIM_LAYER, LEAVES, MEMORY) = spac;

    for _i in 1:DIM_LAYER
        LEAVES[_i].g_H₂O_s_shaded = state.gs_shaded[_i];
        LEAVES[_i].g_H₂O_s_sunlit .= state.gs_sunlit[:,:,_i];
    end;

    MEMORY.tem .= state.t_clm;

    return nothing
);
