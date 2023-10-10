

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-15: generalize the function
#
#######################################################################################################################################################################################################
"""
This function updates the extinction (and scattering) coefficients for canopy. Supported methods are to
- Update coefficients for broadband single layer canopy
- Update coefficients for hyperspectral mutiple layers canopy

"""
function extinction_scattering_coefficients! end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-07: add function to update the extinction and scattering coefficient cache within canopy
#     2022-Jun-07: update _Co, _Cs, _So, _Ss as well for hyperspectral canopy
#     2022-Jun-15: add method for broadband single layer canopy
#     2023-Jun-20: add config to parameter list
#
#######################################################################################################################################################################################################
"""

    extinction_scattering_coefficients!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT}

Update the extinction and scattering coefficients, given
- `config` SPAC configurations
- `can` `MultiLayerCanopy` type canopy

"""
extinction_scattering_coefficients!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT} = (
    (; Θ_INCL) = config;
    (; OPTICS) = can;

    for i in eachindex(Θ_INCL)
        can.sun_geometry.auxil.ks_incl[i],
        OPTICS._ko[i],
        OPTICS._sb[i],
        OPTICS._sf[i],
        OPTICS._Co[i],
        can.sun_geometry.auxil.Cs_incl[i],
        OPTICS._So[i],
        can.sun_geometry.auxil.Ss_incl[i] =
            extinction_coefficient(can.sun_geometry.state.sza, can.sensor_geometry.state.vza, can.sensor_geometry.state.vaa - can.sun_geometry.state.saa, Θ_INCL[i]);
    end;

    return nothing
);
