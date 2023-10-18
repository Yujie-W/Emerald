const MODIS_BAND_1 = [ 620,  670];  # RED
const MODIS_BAND_2 = [ 841,  876];  # NIR
const MODIS_BAND_3 = [ 459,  479];  # BLUE
const MODIS_BAND_7 = [2105, 2155];  # SWIR


#######################################################################################################################################################################################################
#
# Changes to these functions
# General
#     2022-Jun-13: add function to compute MODIS EVI, EVI2, LSWI, NDVI, and NIRv
#     2022-Oct-19: add function to compute MODIS BLUE, NIR, RED, and NIRv radiance
#
#######################################################################################################################################################################################################
"""

    MODIS_BLUE(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return blue band reflectance for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_BLUE end;

MODIS_BLUE(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_BLUE(spac.canopy, config.SPECTRA);

MODIS_BLUE(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ, can.sensor_geometry.auxil.reflectance, FT(MODIS_BAND_3[1]), FT(MODIS_BAND_3[2]); steps=4)
);


"""

    MODIS_EVI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return EVI for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_EVI end;

MODIS_EVI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_EVI(spac.canopy, config.SPECTRA);

MODIS_EVI(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    blue = MODIS_BLUE(can, spectra);
    red  = MODIS_RED(can, spectra);
    nir  = MODIS_NIR(can, spectra);

    return FT(2.5) * (nir - red) / (nir + 6 * red - FT(7.5) * blue + 1)
);


"""

    MODIS_EVI2(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return EVI2 for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_EVI2 end;

MODIS_EVI2(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_EVI2(spac.canopy, config.SPECTRA);

MODIS_EVI2(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    red = MODIS_RED(can, spectra);
    nir = MODIS_NIR(can, spectra);

    return FT(2.5) * (nir - red) / (nir + FT(2.4) * red + 1)
);


"""

    MODIS_LSWI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return LSWI for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_LSWI end;

MODIS_LSWI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_LSWI(spac.canopy, config.SPECTRA);

MODIS_LSWI(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    nir  = MODIS_NIR(can, spectra);
    swir = read_spectrum(spectra.Λ, can.sensor_geometry.auxil.reflectance, FT(MODIS_BAND_7[1]), FT(MODIS_BAND_7[2]); steps=5);

    return (nir - swir) / (nir + swir)
);


"""

    MODIS_NDVI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return NDVI for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_NDVI end;

MODIS_NDVI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_NDVI(spac.canopy, config.SPECTRA);

MODIS_NDVI(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    red = MODIS_RED(can, spectra);
    nir = MODIS_NIR(can, spectra);

    return (nir - red) / (nir + red)
);


"""

    MODIS_NIR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return near infrared band reflectance for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_NIR end;

MODIS_NIR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_NIR(spac.canopy, config.SPECTRA);

MODIS_NIR(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ, can.sensor_geometry.auxil.reflectance, FT(MODIS_BAND_2[1]), FT(MODIS_BAND_2[2]); steps=6)
);


"""

    MODIS_NIRv(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return NIRv for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_NIRv end;

MODIS_NIRv(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_NIRv(spac.canopy, config.SPECTRA);

MODIS_NIRv(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    red = MODIS_RED(can, spectra);
    nir = MODIS_NIR(can, spectra);

    return (nir - red) / (nir + red) * nir
);


"""

    MODIS_NIRvR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return NIRv radiance for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_NIRvR end;

MODIS_NIRvR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_NIRvR(spac.canopy, config.SPECTRA);

MODIS_NIRvR(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    nir_rad = read_spectrum(spectra.Λ, can.sensor_geometry.auxil.e_sensor, FT(MODIS_BAND_2[1]), FT(MODIS_BAND_2[2]); steps=6);

    return MODIS_NDVI(can, spectra) * nir_rad
);


"""

    MODIS_RED(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return red band reflectance for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_RED end;

MODIS_RED(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_RED(spac.canopy, config.SPECTRA);

MODIS_RED(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ, can.sensor_geometry.auxil.reflectance, FT(MODIS_BAND_1[1]), FT(MODIS_BAND_1[2]); steps=6)
);
