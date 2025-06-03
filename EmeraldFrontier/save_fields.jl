#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Aug-06: isolate the function to save the fields to the NamedTuple
#     2024-Aug-06: output the junction pressure
#     2024-Aug-06: read leaf water potential only if the leaf is not shedded; otherwise, set it to NaN
#     2024-Aug-08: save OCS flux if requested
#     2024-Sep-09: save SAP_VOLUME if requested
#
#######################################################################################################################################################################################################
"""

    save_fields!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, wdf::NamedTuple, ind::Int; saving_dict::Dict{String,Any} = SAVING_DICT) where {FT}

Save the fields to the NamedTuple, given
- `config` the configuration of the SPAC model
- `spac` the SPAC model
- `wdf` the NamedTuple to store the outputs
- `ind` the index of the row in the DataFrame
- `saving_dict` the dictionary to store the settings for saving the outputs

"""
function save_fields!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, wdf::NamedTuple, ind::Int; saving_dict::Dict{String,Any} = SAVING_DICT) where {FT}
    # read the data out of dataframe row to reduce memory allocation
    df_dif::FT = wdf.RAD_DIF[ind];
    df_dir::FT = wdf.RAD_DIR[ind];

    # save the profiles of the soil
    if saving_dict["MOD_SWC"]
        for i in eachindex(spac.soils)
            wdf[Symbol("MOD_SWC_$i")][ind] = spac.soils[i].state.θ;
        end;
    end;
    if saving_dict["MOD_P_SOIL"]
        for i in eachindex(spac.soils)
            wdf[Symbol("MOD_P_SOIL_$i")][ind] = spac.soils[i].s_aux.ψ;
        end;
    end;
    if saving_dict["MOD_T_SOIL"]
        for i in eachindex(spac.soils)
            wdf[Symbol("MOD_T_SOIL_$i")][ind] = spac.soils[i].s_aux.t;
        end;
    end;

    # save the profiles of the leaves
    if saving_dict["MOD_T_LEAF"]
        for i in eachindex(spac.plant.leaves)
            wdf[Symbol("MOD_T_LEAF_$i")][ind] = spac.plant.leaves[i].energy.s_aux.t;
        end;
    end;
    if saving_dict["MOD_T_MMM"]
        sum_t::FT = 0;
        min_t::FT = 999;
        max_t::FT = 0;
        for l in spac.plant.leaves
            sum_t += l.energy.s_aux.t;
            min_t = min(min_t, l.energy.s_aux.t);
            max_t = max(max_t, l.energy.s_aux.t);
        end;
        wdf.MOD_T_L_MAX[ind]  = max_t;
        wdf.MOD_T_L_MEAN[ind] = sum_t / length(spac.plant.leaves);
        wdf.MOD_T_L_MIN[ind]  = min_t;
    end;

    # save the CO2 and H2O fluxes
    if saving_dict["BETA"]
        wdf.BETA[ind] = BETA(spac);
    end;
    if saving_dict["CNPP"]
        wdf.CNPP[ind] = CNPP(spac);
    end;
    if saving_dict["ET_SOIL"]
        wdf.ET_SOIL[ind] = ET_SOIL(spac);
    end;
    if saving_dict["ET_VEGE"]
        wdf.ET_VEGE[ind] = ET_VEGE(spac);
    end;
    if saving_dict["GPP"]
        wdf.GPP[ind] = GPP(spac);
    end;
    if saving_dict["OCS"]
        wdf.OCS[ind] = OCS(spac);
    end;

    # save the SIF (PAR and PPAR) if there is sunlight (0 otherwise)
    daytime = df_dir + df_dif >= 10;
    if saving_dict["SIF683"]
        wdf.SIF683[ind] = daytime ? TROPOMI_SIF683(config, spac) : 0;
    end;
    if saving_dict["SIF740"]
        wdf.SIF740[ind] = daytime ? TROPOMI_SIF740(config, spac) : 0;
    end;
    if saving_dict["SIF757"]
        wdf.SIF757[ind] = daytime ? OCO2_SIF759(config, spac) : 0;
    end;
    if saving_dict["SIF771"]
        wdf.SIF771[ind] = daytime ? OCO2_SIF770(config, spac) : 0;
    end;
    if saving_dict["ΣSIF"]
        wdf.ΣSIF[ind] = daytime ? ΣSIF(spac) : 0;
    end;
    if saving_dict["ΣSIF_CHL"]
        wdf.ΣSIF_CHL[ind] = daytime ? ΣSIF_CHL(config, spac) : 0;
    end;
    if saving_dict["ΣSIF_LEAF"]
        wdf.ΣSIF_LEAF[ind] = daytime ? ΣSIF_LEAF(config, spac) : 0;
    end;
    if saving_dict["PAR"]
        wdf.PAR[ind] = daytime ? PAR(config, spac) : 0;
    end;
    if saving_dict["PPAR"]
        wdf.PPAR[ind] = daytime ? PPAR(spac) : 0;
    end;

    # save the VI (and phi) if there is sunlight
    if daytime
        if saving_dict["MOD_ΦDΦN"]
            wdf.ΦF[ind],wdf.ΦP[ind] = ΦD_ΦN(spac);
        end;
        if saving_dict["MOD_ΦFΦP"]
            wdf.ΦF[ind],wdf.ΦP[ind] = ΦF_ΦP(spac);
        end;
        if saving_dict["NDVI"]
            wdf.NDVI[ind] = MODIS_NDVI(config, spac);
        end;
        if saving_dict["EVI"]
            wdf.EVI[ind] = MODIS_EVI(config, spac);
        end;
        if saving_dict["NIRvI"]
            wdf.NIRvI[ind] = MODIS_NIRv(config, spac);
        end;
        if saving_dict["NIRvR"]
            wdf.NIRvR[ind] = MODIS_NIRvR(config, spac);
        end;
    end;

    # save the plant health status
    if saving_dict["C_POOL"]
        wdf.C_POOL[ind] = spac.plant.pool.c_pool;
    end;
    if saving_dict["K_PLANT"]
        wdf.K_PLANT[ind] = K_PLANT(spac);
    end;
    if saving_dict["K_ROOT_STEM"]
        wdf.K_ROOT_STEM[ind] = K_PLANT(spac; include_leaf = false);
    end;
    if saving_dict["MOD_P_LEAF"]
        for i in eachindex(spac.plant.leaves)
            wdf[Symbol("MOD_P_LEAF_$i")][ind] = spac.plant._leaf_shedded ? NaN : spac.plant.leaves[i].xylem.auxil.pressure[end];
        end;
    end;
    if saving_dict["MOD_P_MMM"]
        if spac.plant._leaf_shedded
            wdf.MOD_P_L_MAX[ind]  = NaN;
            wdf.MOD_P_L_MEAN[ind] = NaN;
            wdf.MOD_P_L_MIN[ind]  = NaN;
        else
            sum_p::FT = 0;
            min_p::FT = 0;
            max_p::FT = -999;
            for l in spac.plant.leaves
                sum_p += l.xylem.auxil.pressure[end];
                min_p = min(min_p, l.xylem.auxil.pressure[end]);
                max_p = max(max_p, l.xylem.auxil.pressure[end]);
            end;
            wdf.MOD_P_L_MAX[ind]  = max_p;
            wdf.MOD_P_L_MEAN[ind] = sum_p / length(spac.plant.leaves);
            wdf.MOD_P_L_MIN[ind]  = min_p;
        end;
    end;
    if saving_dict["P_JUNCTION"]
        wdf.P_JUNCTION[ind] = spac.plant.junction.s_aux.pressure;
    end;
    if saving_dict["SAP_VOLUME"]
        wdf.SAP_VOLUME[ind] = SAP_VOLUME(spac);
    end;
    if saving_dict["TRUNK_AREA"]
        wdf.TRUNK_AREA[ind] = spac.plant.trunk.xylem.trait.area;
    end;

    return nothing
end;
