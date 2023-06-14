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
    state.beta = BETA(spac);
    state.csif = ΣSIF(spac);
    state.etr = ΣETR(spac);
    state.gpp = GPP(spac);
    state.modis_evi = MODIS_EVI(CANOPY);
    state.modis_ndvi = MODIS_NDVI(CANOPY);
    state.modis_nirv = MODIS_NIRv(CANOPY);
    state.oco_sif₇₅₉ = OCO2_SIF759(CANOPY);
    state.oco_sif₇₇₀ = OCO2_SIF770(CANOPY);
    state.par = spac.CANOPY.RADIATION.par_in;
    state.ppar = PPAR(spac);
    state.transpiration = T_VEG(spac);
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
