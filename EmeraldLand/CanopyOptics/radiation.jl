#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Apr-13: name the method to shortwave_radiation! to be more specific
#
#######################################################################################################################################################################################################
"""

Run short wave radiation.

"""
function shortwave_radiation! end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-09: migrate the function from CanopyLayers
#     2022-Jun-10: rename PAR/APAR to APAR/PPAR to be more accurate
#     2022-Jun-10: add PAR calculation (before absorption)
#     2022-Jun-10: compute shortwave net radiation
#     2022-Jun-13: use DIM_LAYER instead of _end
#     2022-Jun-29: use Leaf for the hyperspectral RT
#     2023-Mar-11: add code to account for the case of LAI == 0
#     2023-Apr-13: rename option APAR_car to apar_car
#     2023-Apr-13: name the method to shortwave_radiation! to be more specific
#     2023-May-19: use δlai per canopy layer
#     2023-Jun-15: compute PAR when lai = 0
#     2023-Jun-20: remove option apar_car as it is already in config
# Bug fixes
#     2022-Jul-15: sum by r_net_sw by the weights of sunlit and shaded fractions
#     2022-Jul-27: use _ρ_dd, _ρ_sd, _τ_dd, and _τ_sd for leaf energy absorption (typo when refactoring the code)
#     2022-Aug-30: fix par, apar, and ppar issues
#     2022-Nov-29: ppar for sunlit leaves should be scaled based on sunlit fraction
#
#######################################################################################################################################################################################################
"""

    shortwave_radiation!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::ShortwaveRadiation{FT}, sbulk::SoilBulk{FT}) where {FT}

Updates canopy radiation profiles for shortwave radiation, given
- `config` Configuration for `MultiLayerSPAC`
- `can` `MultiLayerCanopy` type struct
- `leaves` Vector of `Leaf`
- `rad` Incoming shortwave radiation
- `sbulk` Soil bulk parameters

"""
shortwave_radiation!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::ShortwaveRadiation{FT}, sbulk::SoilBulk{FT}) where {FT} = (
    (; DIM_LAYER) = config;
    (; OPTICS, RADIATION) = can;

    # 3. compute the spectra at the observer direction
    for i in 1:DIM_LAYER
        _e_d_i = view(RADIATION.e_diffuse_down,:,i);   # downward diffuse radiation at upper boundary
        _e_u_i = view(RADIATION.e_diffuse_up  ,:,i);   # upward diffuse radiation at upper boundary
        _e_v_i = view(RADIATION.e_v,:,i);

        _dob_i = view(OPTICS.σ_dob,:,i);   # scattering coefficient backward for diffuse->observer
        _dof_i = view(OPTICS.σ_dof,:,i);   # scattering coefficient forward for diffuse->observer
        _so__i = view(OPTICS.σ_so ,:,i);   # bidirectional from solar to observer

        _e_v_i  .= can.sensor_geometry.auxil.po[i] .* _dob_i .* _e_d_i;
        _e_v_i .+= can.sensor_geometry.auxil.po[i] .* _dof_i .* _e_u_i;
        _e_v_i .+= can.sensor_geometry.auxil.pso[i] .* _so__i .* rad.e_dir;
        _e_v_i .*= can.structure.state.δlai[i] * can.structure.auxil.ci;
    end;
    RADIATION.e_v[:,end] .= can.sensor_geometry.auxil.po[end] .* view(RADIATION.e_diffuse_up,:,DIM_LAYER+1);

    for i in eachindex(RADIATION.e_o)
        RADIATION.e_o[i] = sum(view(RADIATION.e_v,i,:)) / FT(pi);
    end;

    RADIATION.albedo .= RADIATION.e_o * FT(pi) ./ (rad.e_dir .+ rad.e_dif);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Apr-13: move method of canopy_radiation! to longwave_radiation!
#
#######################################################################################################################################################################################################
"""

Update long wave radiation profiles for canopy.

"""
function longwave_radiation! end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-10: migrate the function thermal_fluxes! from CanopyLayers
#     2022-Jun-10: update net lw radiation for leaves and soil
#     2022-Jun-29: use Leaf for the hyperspectral RT
#     2023-Mar-11: add code to account for the case of LAI == 0
#
#######################################################################################################################################################################################################
"""

    longwave_radiation!(can::MultiLayerCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::FT, sbulk::SoilBulk{FT}, soil::SoilLayer{FT}) where {FT}

Updates canopy radiation profiles for shortwave or longwave radiation, given
- `can` `MultiLayerCanopy` type struct
- `leaves` Vector of `Leaf`
- `rad` Incoming longwave radiation
- `sbulk` Soil bulk parameters
- `soil` First soil layer

"""
longwave_radiation!(can::MultiLayerCanopy{FT}, leaves::Vector{Leaf{FT}}, rad::FT, sbulk::SoilBulk{FT}, soil::SoilLayer{FT}) where {FT} = (

    if can.structure.state.lai == 0
        _r_lw_soil = K_STEFAN(FT) * (1 - sbulk.auxil.ρ_lw) * soil.auxil.t ^ 4;
        RADIATION.r_lw .= 0;
        RADIATION.r_net_lw .= 0;
        RADIATION.r_lw_up .= rad * sbulk.auxil.ρ_lw + _r_lw_soil;
        sbulk.auxil.r_net_lw = rad * (1 - sbulk.auxil.ρ_lw) - _r_lw_soil;

        return nothing
    end;

    return nothing
);


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
