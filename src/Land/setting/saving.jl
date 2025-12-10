""" Default saving settings """
DEFAULT_SAVING_SETTINGS = ParameterFunctionMapper[
    ParameterFunctionMapper("APAR", false, APAR, []),
    ParameterFunctionMapper("PAR", false, PAR, []),
    ParameterFunctionMapper("PPAR", false, PPAR, []),

    ParameterFunctionMapper("BETA", false, BETA, []),
    ParameterFunctionMapper("PCI", false, LEAF_PCI, []),

    ParameterFunctionMapper("CNPP", false, CNPP, []),
    ParameterFunctionMapper("ET", true, ET, []),
    ParameterFunctionMapper("ET_SOIL", false, ET_SOIL, []),
    ParameterFunctionMapper("ET_VEGE", false, ET_VEGE, []),
    ParameterFunctionMapper("GPP", true, GPP, []),
    ParameterFunctionMapper("OCS", false, OCS, []),

    ParameterFunctionMapper("EVI", false, MODIS_EVI, []),
    ParameterFunctionMapper("EVI2", false, MODIS_EVI2, []),
    ParameterFunctionMapper("LSWI", false, MODIS_LSWI, []),
    ParameterFunctionMapper("NDVI", false, MODIS_NDVI, []),
    ParameterFunctionMapper("NIRv", false, MODIS_NIRv, []),
    ParameterFunctionMapper("NIRvR", false, MODIS_NIRvR, []),

    ParameterFunctionMapper("LATENT_HEAT", false, LATENT_HEAT, []),
    ParameterFunctionMapper("LONGWAVE_OUT", false, LONGWAVE_OUT, []),
    ParameterFunctionMapper("NET_LONGWAVE", false, NET_LONGWAVE, []),
    ParameterFunctionMapper("NET_SHORTWAVE", false, NET_SHORTWAVE, []),
    ParameterFunctionMapper("SENSIBLE_HEAT", false, SENSIBLE_HEAT, []),
    ParameterFunctionMapper("SHORTWAVE_OUT", false, SHORTWAVE_OUT, []),
    ParameterFunctionMapper("T_SKIN", false, T_SKIN, []),

    ParameterFunctionMapper("SIF683", false, TROPOMI_SIF683, []),
    ParameterFunctionMapper("SIF740", true, TROPOMI_SIF740, []),
    ParameterFunctionMapper("SIF747", false, TROPOMI_SIF747, []),
    ParameterFunctionMapper("SIF757", false, OCO2_SIF759, []),
    ParameterFunctionMapper("SIF771", false, OCO2_SIF770, []),
    ParameterFunctionMapper("ΣSIF", false, ΣSIF, []),
    ParameterFunctionMapper("ΣSIF_CHL", false, ΣSIF_CHL, []),
    ParameterFunctionMapper("ΣSIF_LEAF", false, ΣSIF_LEAF, []),

    ParameterFunctionMapper("ΦD", false, ΦD, []),
    ParameterFunctionMapper("ΦF", false, ΦF, []),
    ParameterFunctionMapper("ΦN", false, ΦN, []),
    ParameterFunctionMapper("ΦP", false, ΦP, []),
];


"""

    parameters_to_save(varnames::Vector{String} = String[]; save_all::Bool = false)

Create a saving setting vector for simulation, given
- `varnames` Vector of variable names to be saved
- `save_all` If true, set all parameters to be saved

"""
function parameters_to_save(varnames::Vector{String} = String[]; save_all::Bool = false)
    new_vec = deepcopy(DEFAULT_SAVING_SETTINGS);

    # If save_all is true, set all values to true
    if save_all
        for lpm in new_vec
            lpm.to_save = true;
        end;

        return new_vec
    end;

    # otherwise, loop through the varnames and set the corresponding keys to true
    for vn in varnames
        for lpm in new_vec
            if lpm.name == vn
                lpm.to_save = true;
                break;
            end;
        end;
    end;

    return new_vec
end;
