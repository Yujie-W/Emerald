# All keys start with "MOD_" results in multiple outputs, deal with them carefully in the prepare_df.jl file
const DEFAULT_SAVING_DICT = Dict{String,Bool}(
    # Modeled soil water content and temperature
            "MOD_SWC"     => true,
            "MOD_P_SOIL"  => true,
            "MOD_T_SOIL"  => true,
    # Modeled leaf temperature
            "MOD_T_LEAF"  => false,
            "MOD_T_MMM"   => true,
    # Modeled CO2, H2O, and OCS fluxes
            "BETA"        => false,
            "CNPP"        => true,
            "ET_SOIL"     => true,
            "ET_VEGE"     => true,
            "GPP"         => true,
            "OCS"         => true,
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
            "C_POOL"      => true,
            "K_PLANT"     => true,
            "K_ROOT_STEM" => true,
            "MOD_P_LEAF"  => false,
            "MOD_P_MMM"   => true,
            "P_JUNCTION"  => true,
            "SAP_VOLUME"  => true,
            "TRUNK_AREA"  => true,
    # Modeled heat fluxes
            "MOD_HEAT"    => true,
);


"""

    parameters_to_save(; save_all::Bool = false)

Create a saving dict for simulation, given
- `save_all` If true, set all parameters to be saved

"""
function parameters_to_save(; save_all::Bool = false)
    new_dict = deepcopy(DEFAULT_SAVING_DICT);

    # If save_all is true, set all values to true
    if save_all
        for (k, _) in new_dict
            new_dict[k] = true;
        end;
    end;

    return new_dict
end;
