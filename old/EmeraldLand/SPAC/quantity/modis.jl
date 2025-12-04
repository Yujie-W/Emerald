const MODIS_BAND_1 = [ 620,  670];  # RED
const MODIS_BAND_2 = [ 841,  876];  # NIR
const MODIS_BAND_3 = [ 459,  479];  # BLUE
const MODIS_BAND_4 = [ 545,  565];  # GREEN
const MODIS_BAND_5 = [1230, 1250];  # SWIR
const MODIS_BAND_6 = [1628, 1652];  # SWIR
const MODIS_BAND_7 = [2105, 2155];  # SWIR
const MODIS_BANDS = [MODIS_BAND_1, MODIS_BAND_2, MODIS_BAND_3, MODIS_BAND_4, MODIS_BAND_5, MODIS_BAND_6, MODIS_BAND_7];


#######################################################################################################################################################################################################
#
# Changes to these functions
# General
#     2024-Jul-10: add function to compute MODIS band reflectance (general)
#     2024-Jul-10: add option to weight the reflectance based on radiation
#
#######################################################################################################################################################################################################
"""

    MODIS_BAND_REFL(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, i::Int; steps::Int = 4) where {FT}

Return band reflectance for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC
- `i` Band index
- `steps` Number of steps to compute reflectance
- `weighted` Weigh reflectance based on radiation

"""
function MODIS_BAND_REFL end;

MODIS_BAND_REFL(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, i::Int; steps::Int = 4, weighted::Bool = false) where {FT} =
    MODIS_BAND_REFL(spac.canopy, config.SPECTRA, spac.meteo.rad_sw, i; steps = steps, weighted = weighted);

MODIS_BAND_REFL(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}, radiation::ShortwaveRadiation{FT}, i::Int; steps::Int = 4, weighted::Bool = false) where {FT} = (
    if weighted
        return read_spectrum(spectra.Λ, can.sensor_geometry.auxil.reflectance, radiation.e_dir .+ radiation.e_dif, FT(MODIS_BANDS[i][1]), FT(MODIS_BANDS[i][2]); steps = steps)
    end;

    return read_spectrum(spectra.Λ, can.sensor_geometry.auxil.reflectance, FT(MODIS_BANDS[i][1]), FT(MODIS_BANDS[i][2]); steps = steps)
);

MODIS_RED(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT} = MODIS_BAND_REFL(config, spac, 1; steps = 6, weighted = weighted);
MODIS_NIR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT} = MODIS_BAND_REFL(config, spac, 2; steps = 6, weighted = weighted);
MODIS_BLUE(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT} = MODIS_BAND_REFL(config, spac, 3; steps = 4, weighted = weighted);
MODIS_SWIR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT} = MODIS_BAND_REFL(config, spac, 7; steps = 5, weighted = weighted);


#######################################################################################################################################################################################################
#
# Changes to these functions
# General
#     2022-Jun-13: add function to compute MODIS EVI, EVI2, LSWI, NDVI, and NIRv
#     2022-Oct-19: add function to compute MODIS BLUE, NIR, RED, and NIRv radiance
#     2024-Jul-10: add option to weight the reflectance based on radiation
#
#######################################################################################################################################################################################################
"""

    MODIS_EVI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}

Return EVI for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC
- `weighted` Weigh reflectance based on radiation

"""
function MODIS_EVI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}
    blue = MODIS_BLUE(config, spac; weighted = weighted);
    red  = MODIS_RED(config, spac; weighted = weighted);
    nir  = MODIS_NIR(config, spac; weighted = weighted);

    return FT(2.5) * (nir - red) / (nir + 6 * red - FT(7.5) * blue + 1)
end;


"""

    MODIS_EVI2(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}

Return EVI2 for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC
- `weighted` Weigh reflectance based on radiation

"""
function MODIS_EVI2(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}
    red = MODIS_RED(config, spac; weighted = weighted);
    nir = MODIS_NIR(config, spac; weighted = weighted);

    return FT(2.5) * (nir - red) / (nir + FT(2.4) * red + 1)
end;


"""

    MODIS_LSWI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}

Return LSWI for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC
- `weighted` Weigh reflectance based on radiation

"""
function MODIS_LSWI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}
    nir  = MODIS_NIR(config, spac; weighted = weighted);
    swir = MODIS_SWIR(config, spac; weighted = weighted);

    return (nir - swir) / (nir + swir)
end;


"""

    MODIS_NDVI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}

Return NDVI for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC
- `weighted` Weigh reflectance based on radiation

"""
function MODIS_NDVI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}
    red = MODIS_RED(config, spac; weighted = weighted);
    nir = MODIS_NIR(config, spac; weighted = weighted);

    return (nir - red) / (nir + red)
end;


"""

    MODIS_NIRv(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}

Return NIRv for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC
- `weighted` Weigh reflectance based on radiation

"""
function MODIS_NIRv(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}
    red = MODIS_RED(config, spac; weighted = weighted);
    nir = MODIS_NIR(config, spac; weighted = weighted);

    return (nir - red) / (nir + red) * nir
end;


"""

    MODIS_NIRvR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}

Return NIRv radiance for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC
- `weighted` Weigh reflectance based on radiation

"""
function MODIS_NIRvR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}; weighted::Bool = false) where {FT}
    nir_rad = read_spectrum(config.SPECTRA.Λ, spac.canopy.sensor_geometry.auxil.e_sensor, FT(MODIS_BAND_2[1]), FT(MODIS_BAND_2[2]); steps=6);

    return MODIS_NDVI(config, spac; weighted = weighted) * nir_rad
end;
