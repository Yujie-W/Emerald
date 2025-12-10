"""

    save_fields!(config::SPACConfig{FT}, spac::BulkSPAC{FT}, results::NamedTuple, ind::Int; saving_setting::Dict{String,Bool} = SAVING_DICT) where {FT}

Save the fields to the NamedTuple, given
- `config` the configuration of the SPAC model
- `spac` the SPAC model
- `results` the NamedTuple to store the outputs
- `ind` the index of the row in the DataFrame
- `saving_setting` the dictionary to store the settings for saving the outputs

"""
function save_fields!(config::SPACConfig{FT}, spac::BulkSPAC{FT}, results::NamedTuple, ind::Int; saving_setting::Vector{ParameterFunctionMapper} = parameters_to_save()) where {FT}
    for lpm in saving_setting
        if lpm.to_save
            results[Symbol(lpm.name)][ind] = lpm.func(config, spac, lpm.params...);
        end;
    end;
    #=
    # save the profiles of the soil
    if saving_setting["MOD_SWC"]
        for i in eachindex(spac.soils)
            results[Symbol("MOD_SWC_$i")][ind] = spac.soils[i].state.θ;
            results[Symbol("MOD_SWC_ICE_$i")][ind] = spac.soils[i].state.θ_ice;
        end;
    end;
    if saving_setting["MOD_P_SOIL"]
        for i in eachindex(spac.soils)
            results[Symbol("MOD_P_SOIL_$i")][ind] = spac.soils[i].s_aux.ψ;
        end;
    end;
    if saving_setting["MOD_T_SOIL"]
        for i in eachindex(spac.soils)
            results[Symbol("MOD_T_SOIL_$i")][ind] = spac.soils[i].s_aux.t;
        end;
    end;

    # save the profiles of the leaves
    if saving_setting["MOD_T_LEAF"]
        for i in eachindex(spac.plant.leaves)
            results[Symbol("MOD_T_LEAF_$i")][ind] = spac.plant.leaves[i].energy.s_aux.t;
        end;
    end;
    if saving_setting["MOD_T_MMM"]
        sum_t::FT = 0;
        min_t::FT = 999;
        max_t::FT = 0;
        for l in spac.plant.leaves
            sum_t += l.energy.s_aux.t;
            min_t = min(min_t, l.energy.s_aux.t);
            max_t = max(max_t, l.energy.s_aux.t);
        end;
        results.MOD_T_L_MAX[ind]  = max_t;
        results.MOD_T_L_MEAN[ind] = sum_t / length(spac.plant.leaves);
        results.MOD_T_L_MIN[ind]  = min_t;
    end;

    # save the CO2 and H2O fluxes
    if saving_setting["CNPP"]
        results.CNPP[ind] = CNPP(spac);
    end;
    if saving_setting["ET_SOIL"]
        results.ET_SOIL[ind] = ET_SOIL(spac);
    end;
    if saving_setting["ET_VEGE"]
        results.ET_VEGE[ind] = ET_VEGE(spac);
    end;
    if saving_setting["PCI"]
        results.PCI[ind] = LEAF_PCI(spac);
    end;

    # save the SIF (PAR and PPAR) if there is sunlight (0 otherwise)
    daytime = PAR(config, spac) > 0;
    if saving_setting["SIF683"]
        results.SIF683[ind] = daytime ? TROPOMI_SIF683(config, spac) : 0;
    end;
    if saving_setting["SIF740"]
        results.SIF740[ind] = daytime ? TROPOMI_SIF740(config, spac) : 0;
    end;
    if saving_setting["SIF757"]
        results.SIF757[ind] = daytime ? OCO2_SIF759(config, spac) : 0;
    end;
    if saving_setting["SIF771"]
        results.SIF771[ind] = daytime ? OCO2_SIF770(config, spac) : 0;
    end;
    if saving_setting["ΣSIF"]
        results.ΣSIF[ind] = daytime ? ΣSIF(config, spac) : 0;
    end;
    if saving_setting["ΣSIF_CHL"]
        results.ΣSIF_CHL[ind] = daytime ? ΣSIF_CHL(config, spac) : 0;
    end;
    if saving_setting["ΣSIF_LEAF"]
        results.ΣSIF_LEAF[ind] = daytime ? ΣSIF_LEAF(config, spac) : 0;
    end;
    if saving_setting["PAR"]
        results.PAR[ind] = daytime ? PAR(config, spac) : 0;
    end;
    if saving_setting["APAR"]
        results.APAR[ind] = daytime ? APAR(spac) : 0;
    end;
    if saving_setting["PPAR"]
        results.PPAR[ind] = daytime ? PPAR(spac) : 0;
    end;

    # save the VI (and phi) if there is sunlight
    if daytime
        if saving_setting["MOD_ΦDΦN"]
            results.ΦF[ind],results.ΦP[ind] = ΦD_ΦN(spac);
        end;
        if saving_setting["MOD_ΦFΦP"]
            results.ΦF[ind],results.ΦP[ind] = ΦF_ΦP(spac);
        end;
        if saving_setting["NDVI"]
            results.NDVI[ind] = MODIS_NDVI(config, spac);
        end;
        if saving_setting["EVI"]
            results.EVI[ind] = MODIS_EVI(config, spac);
        end;
        if saving_setting["NIRvI"]
            results.NIRvI[ind] = MODIS_NIRv(config, spac);
        end;
        if saving_setting["NIRvR"]
            results.NIRvR[ind] = MODIS_NIRvR(config, spac);
        end;
    end;

    # save the plant health status
    if saving_setting["C_POOL"]
        results.C_POOL[ind] = spac.plant.pool.c_pool;
    end;
    if saving_setting["K_PLANT"]
        results.K_PLANT[ind] = K_PLANT(spac);
    end;
    if saving_setting["K_ROOT_STEM"]
        results.K_ROOT_STEM[ind] = K_PLANT(spac; include_leaf = false);
    end;
    if saving_setting["MOD_P_LEAF"]
        for i in eachindex(spac.plant.leaves)
            results[Symbol("MOD_P_LEAF_$i")][ind] = spac.plant._leaf_shedded ? NaN : spac.plant.leaves[i].xylem.auxil.pressure[end];
        end;
    end;
    if saving_setting["MOD_P_MMM"]
        if spac.plant._leaf_shedded
            results.MOD_P_L_MAX[ind]  = NaN;
            results.MOD_P_L_MEAN[ind] = NaN;
            results.MOD_P_L_MIN[ind]  = NaN;
        else
            sum_p::FT = 0;
            min_p::FT = 0;
            max_p::FT = -999;
            for l in spac.plant.leaves
                sum_p += l.xylem.auxil.pressure[end];
                min_p = min(min_p, l.xylem.auxil.pressure[end]);
                max_p = max(max_p, l.xylem.auxil.pressure[end]);
            end;
            results.MOD_P_L_MAX[ind]  = max_p;
            results.MOD_P_L_MEAN[ind] = sum_p / length(spac.plant.leaves);
            results.MOD_P_L_MIN[ind]  = min_p;
        end;
    end;
    if saving_setting["P_JUNCTION"]
        results.P_JUNCTION[ind] = spac.plant.junction.s_aux.pressure;
    end;
    if saving_setting["SAP_VOLUME"]
        results.SAP_VOLUME[ind] = SAP_VOLUME(spac);
    end;
    if saving_setting["TRUNK_AREA"]
        results.TRUNK_AREA[ind] = spac.plant.trunk.xylem.trait.area;
    end;

    # save the heat fluxes
    if saving_setting["MOD_HEAT"]
        results.MOD_LATENT_HEAT[ind] = LATENT_HEAT(spac);
        results.MOD_SENSIBLE_HEAT[ind] = SENSIBLE_HEAT(spac);
        results.MOD_NET_LONGWAVE[ind] = NET_LONGWAVE(spac);
        results.MOD_NET_SHORTWAVE[ind] = NET_SHORTWAVE(spac);
        results.MOD_LONGWAVE_OUT[ind] = LONGWAVE_OUT(spac);
        results.MOD_SHORTWAVE_OUT[ind] = SHORTWAVE_OUT(config, spac);
    end;
    =#

    return nothing
end;
