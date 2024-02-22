#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to synchronize state among spac and spac state structs
#
#######################################################################################################################################################################################################
"""

    spac_state!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, state::MultiLayerSPACState{FT}) where {FT}
    spac_state!(state::MultiLayerSPACState{FT}, spac::BulkSPAC{FT}) where {FT}

Synchronize state variables from 1st to 2nd struct, given
- `config` SPAC configurations
- `spac` `BulkSPAC` struct for SPAC
- `state` `MultiLayerSPACState` struct for states

"""
function spac_state! end;

spac_state!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, state::MultiLayerSPACState{FT}) where {FT} = (
    leaves = spac.plant.leaves;
    memory = spac.plant.memory;

    for i in eachindex(leaves)
        state.gs_shaded[i] = leaves[i].g_H₂O_s_shaded;
        state.gs_sunlit[:,:,i] .= leaves[i].g_H₂O_s_sunlit;
    end;

    state.t_clm .= memory.t_history;

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
    state.par = PAR(config, spac);
    state.ppar = PPAR(spac);
    state.transpiration = T_VEG(spac);
    state.tropomi_sif₆₈₃ = TROPOMI_SIF683(config, spac);
    state.tropomi_sif₇₄₀ = TROPOMI_SIF740(config, spac);

    return nothing
);

spac_state!(state::MultiLayerSPACState{FT}, spac::BulkSPAC{FT}) where {FT} = (
    leaves = spac.plant.leaves;
    memory = spac.plant.memory;

    for i in eachindex(leaves)
        leaves[i].g_H₂O_s_shaded = state.gs_shaded[i];
        leaves[i].g_H₂O_s_sunlit .= state.gs_sunlit[:,:,i];
    end;

    memory.t_history .= state.t_clm;

    return nothing
);
