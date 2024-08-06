# Netcdf configuration for output
SAVING_DICT = Dict{String, Any}(
    # Modeled soil water content and temperature
            "MOD_SWC"     => true,
            "MOD_P_SOIL"  => true,
            "MOD_T_SOIL"  => true,
    # Modeled leaf temperature
            "MOD_T_LEAF"  => false,
            "MOD_T_MMM"   => true,
    # Modeled CO2 and H2O fluxes
            "BETA"        => false,
            "CNPP"        => true,
            "GPP"         => true,
            "ET_VEG"      => true,
    # SIF (default is false)
            "SIF683"      => false,
            "SIF740"      => true,
            "SIF757"      => false,
            "SIF771"      => false,
            "ΣSIF"        => false,
            "ΣSIF_CHL"    => false,
            "ΣSIF_LEAF"   => false,
            "ΦFΦP"        => false,
    # VI (default is false)
            "NDVI"        => false,
            "EVI"         => false,
            "NIRvI"       => false,
            "NIRvR"       => false,
            "PAR"         => false,
            "PPAR"        => false,
    # Modeled plant health status
            "K_PLANT"     => true,
            "K_ROOT_STEM" => true,
            "MOD_P_LEAF"  => false,
            "MOD_P_MMM"   => true,
            "P_JUNCTION"  => true,
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Apr-13: add function to create spac configuration
#
#######################################################################################################################################################################################################
"""

    spac_config(gm_dict::Dict)

Create a SPAC configuration struct, given
- `gm_dict` Dictionary of GriddingMachine data in a grid

"""
function spac_config(gm_dict::Dict)
    config = SPACConfiguration(gm_dict["FT"]);
    config.MESSAGE_LEVEL = gm_dict["MESSAGE_LEVEL"];
    config.DIM_PPAR_BINS = 10;

    return config
end;
