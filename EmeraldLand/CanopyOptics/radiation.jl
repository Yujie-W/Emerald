




#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-09: migrate the function from CanopyLayers
#     2022-Jun-09: rename function to canopy_radiation!
#
#######################################################################################################################################################################################################
"""

Run shortwave and longwave radiation.

"""
function canopy_radiation! end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-29: add method for SPAC
#     2022-Jul-28: update soil albedo at the very first step
#     2023-Mar-11: run canopy optical properties and shortwave radiation only if solar zenith angle is lower than 89
#     2023-Apr-13: sw and lw radiation moved to METEO
#     2023-Jun-15: set albedo to NaN when sza >= 90
#
#######################################################################################################################################################################################################
"""

    canopy_radiation!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Updates canopy radiation profiles for shortwave and longwave radiation, given
- `config` Configurations of spac model
- `spac` `MultiLayerSPAC` type SPAC

"""
canopy_radiation!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; CANOPY, LEAVES, METEO, SOIL_BULK, SOILS) = spac;
    (; DIM_LAYER) = config;

    soil_albedo!(config, SOIL_BULK, SOILS[1]);
    # TODO: note here that this will disable the optical properties of longwave radiation and result in bugs
    if CANOPY.sun_geometry.state.sza < 89
        canopy_optical_properties!(config, CANOPY);
        canopy_optical_properties!(config, CANOPY, LEAVES, SOIL_BULK);
        shortwave_radiation!(config, CANOPY, LEAVES, METEO.rad_sw, SOIL_BULK);
    else
        CANOPY.RADIATION.r_net_sw .= 0;
        SOIL_BULK.auxil.r_net_sw = 0;
        CANOPY.RADIATION.par_in_diffuse = 0;
        CANOPY.RADIATION.par_in_direct = 0;
        CANOPY.RADIATION.par_in = 0;
        CANOPY.RADIATION.par_shaded .= 0;
        CANOPY.RADIATION.par_sunlit .= 0;
        CANOPY.RADIATION.apar_shaded .= 0;
        CANOPY.RADIATION.apar_sunlit .= 0;
        CANOPY.RADIATION.e_v .= 0;
        CANOPY.RADIATION.e_o .= 0;
        CANOPY.RADIATION.albedo .= NaN;

        for i in 1:DIM_LAYER
            # PPAR for leaves
            LEAVES[i].flux.auxil.ppar_shaded = 0;
            LEAVES[i].flux.auxil.ppar_sunlit .= 0;
        end;
    end;
    longwave_radiation!(CANOPY, LEAVES, METEO.rad_lw, SOIL_BULK, SOILS[1]);

    return nothing
);
