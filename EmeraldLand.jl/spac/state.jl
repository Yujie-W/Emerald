#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to synchronize state among spac and spac state structs
#
#######################################################################################################################################################################################################
"""

    spac_state!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, state::MultiLayerSPACState{FT}) where {FT}
    spac_state!(state::MultiLayerSPACState{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Synchronize state variables from 1st to 2nd struct, given
- `config` SPAC configurations
- `spac` `MultiLayerSPAC` struct for SPAC
- `state` `MultiLayerSPACState` struct for states

"""
function spac_state! end

spac_state!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, state::MultiLayerSPACState{FT}) where {FT} = (
    (; LEAVES, MEMORY) = spac;

    for _i in eachindex(LEAVES)
        state.gs_shaded[_i] = LEAVES[_i].g_H₂O_s_shaded;
        state.gs_sunlit[:,:,_i] .= LEAVES[_i].g_H₂O_s_sunlit;
    end;

    state.t_clm .= MEMORY.tem;

    # save the variables used for publications
    state.beta = BETA(spac);
    state.csif = ΣSIF(spac);
    state.etr = ΣETR(spac);
    state.gpp = GPP(spac);
    state.modis_evi = MODIS_EVI(config, spac);
    state.modis_ndvi = MODIS_NDVI(config, spac);
    state.modis_nirv = MODIS_NIRv(config, spac);
    state.oco_sif₇₅₉ = OCO2_SIF759(config, spac);
    state.oco_sif₇₇₀ = OCO2_SIF770(config, spac);
    state.par = spac.CANOPY.RADIATION.par_in;
    state.ppar = PPAR(spac);
    state.transpiration = T_VEG(spac);
    state.tropomi_sif₆₈₃ = TROPOMI_SIF683(config, spac);
    state.tropomi_sif₇₄₀ = TROPOMI_SIF740(config, spac);

    return nothing
);

spac_state!(state::MultiLayerSPACState{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; LEAVES, MEMORY) = spac;

    for _i in eachindex(LEAVES)
        LEAVES[_i].g_H₂O_s_shaded = state.gs_shaded[_i];
        LEAVES[_i].g_H₂O_s_sunlit .= state.gs_sunlit[:,:,_i];
    end;

    MEMORY.tem .= state.t_clm;

    return nothing
);
