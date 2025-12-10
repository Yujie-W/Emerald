""" Default saving settings """
DEFAULT_SAVING_SETTINGS = ParameterFunctionMapper[
    ParameterFunctionMapper("APAR", false, APAR, []),
    ParameterFunctionMapper("BETA", false, BETA, []),
    ParameterFunctionMapper("CNPP", false, CNPP, []),
    ParameterFunctionMapper("ET", true, ET, []),
    ParameterFunctionMapper("ET_SOIL", false, ET_SOIL, []),
    ParameterFunctionMapper("ET_VEGE", false, ET_VEGE, []),
    ParameterFunctionMapper("GPP", true, GPP, []),
    ParameterFunctionMapper("OCS", false, OCS, []),
    ParameterFunctionMapper("PCI", false, LEAF_PCI, []),
    ParameterFunctionMapper("PAR", false, PAR, []),
    ParameterFunctionMapper("PPAR", false, PPAR, []),
];


#=
const DEFAULT_SAVING_DICT = Dict{String,Bool}(
    # Modeled soil water content and temperature
            "MOD_SWC"     => true,
            "MOD_P_SOIL"  => false,
            "MOD_T_SOIL"  => true,
    # Modeled leaf temperature
            "MOD_T_LEAF"  => false,
            "MOD_T_MMM"   => false,
    # SIF (default is false)
            "SIF683"      => false,
            "SIF740"      => true,
            "SIF757"      => false,
            "SIF771"      => false,
            "ΣSIF"        => false,
            "ΣSIF_CHL"    => false,
            "ΣSIF_LEAF"   => false,
            "MOD_ΦDΦN"    => false,
            "MOD_ΦFΦP"    => false,
    # VI (default is false)
            "NDVI"        => false,
            "EVI"         => false,
            "NIRvI"       => false,
            "NIRvR"       => false,
            "PAR"         => false,
            "APAR"        => false,
            "PPAR"        => false,
    # Modeled plant health status
            "C_POOL"      => false,
            "K_PLANT"     => false,
            "K_ROOT_STEM" => false,
            "MOD_P_LEAF"  => false,
            "MOD_P_MMM"   => false,
            "P_JUNCTION"  => false,
            "SAP_VOLUME"  => false,
            "TRUNK_AREA"  => false,
    # Modeled heat fluxes
            "MOD_HEAT"    => true,
);
=#


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
