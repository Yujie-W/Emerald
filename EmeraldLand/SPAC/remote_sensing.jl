# Satellite wavelength windows
const MODIS_BAND_1 = [ 620,  670];  # RED
const MODIS_BAND_2 = [ 841,  876];  # NIR
const MODIS_BAND_3 = [ 459,  479];  # BLUE
const MODIS_BAND_7 = [2105, 2155];  # SWIR

const OCO2_SIF_759 = [758.17, 759.20];  # SIF 757
const OCO2_SIF_770 = [769.62, 770.28];  # SIF 770
const OCO3_SIF_759 = [758.26, 759.28];  # SIF 757
const OCO3_SIF_770 = [769.67, 770.34];  # SIF 770

const TROPOMI_SIF_683 = [680, 685];     # SIF 683
const TROPOMI_SIF_747 = [735, 758];     # SIF 747
const TROPOMI_SIF_751 = [743, 758];     # SIF 751


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-13: add function to interpolate the spectrum
#
#######################################################################################################################################################################################################
"""
This function interpolate the spectrum to give values at the target wavelength bin(s). The supported methods include
- Interpolate the spectrum at a given wavelength
- Interpolate the spectrum in a given wavelength range

"""
function read_spectrum end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-13: add method to interpolate the spectrum
#
#######################################################################################################################################################################################################
"""

    read_spectrum(x::Vector{FT}, y::Vector{FT}, target::FT) where {FT}

Return the spectrum value at target wavelength bin, given
- `x` X-axis of the spectrum
- `y` Y-axis of the spectrum
- `target` Target x value

"""
read_spectrum(x::Vector{FT}, y::Vector{FT}, target::FT) where {FT} = (
    @assert length(x) == length(y) "Dimensions of provided spectrum x and y must match!";
    @assert x[1] <= target <= x[end] "Target wavelength must be within the range provided spectum!";

    # iterate through the spectrum and find the index
    _ind = 0;
    for i in 1:length(x)-1
        if x[i] <= target <= x[i+1]
            _ind = i;
            break;
        end;
    end;

    return ((x[_ind+1] - target) * y[_ind] + (target - x[_ind]) * y[_ind+1]) / (x[_ind+1] - x[_ind])
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-13: add method to interpolate the spectrum via multiple steps
#
#######################################################################################################################################################################################################
"""

    read_spectrum(x::Vector{FT}, y::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT}

Return the spectrum value at target wavelength bin, given
- `x` X-axis of the spectrum
- `y` Y-axis of the spectrum
- `x₁` Lower x boundary
- `x₂` Upper x boundary
- `steps` The incremental Δx is `(x₂ - x₁) / steps`

"""
read_spectrum(x::Vector{FT}, y::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT} = (
    _ys = 0;
    _δx = (x₂ - x₁) / steps;
    for i in 1:(steps+1)
        _x = x₁ + (i - 1) * _δx;
        _ys += read_spectrum(x, y, _x);
    end;

    return _ys / (steps + 1)
);


#######################################################################################################################################################################################################
#
# Changes to these functions
# General
#     2022-Jun-13: add function to compute MODIS EVI, EVI2, LSWI, NDVI, and NIRv
#     2022-Jun-13: add function to compute OCO2 SIF @ 758.7 and 770.0 nm
#     2022-Jun-13: add function to compute OCO3 SIF @ 758.8 and 770.0 nm
#     2022-Jun-13: add function to compute TROPOMI SIF @ 682.5, 740.0, 746.5, and 750.5 nm
#     2022-Oct-19: add function to compute MODIS BLUE, NIR, RED, and NIRv radiance
#     2023-Jun-16: add methods to compute the quantities directly from spac
#
#######################################################################################################################################################################################################
"""

    MODIS_BLUE(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return blue band reflectance for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_BLUE end;

MODIS_BLUE(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_BLUE(spac.CANOPY, config.SPECTRA);

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

MODIS_EVI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_EVI(spac.CANOPY, config.SPECTRA);

MODIS_EVI(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    _blue = MODIS_BLUE(can, spectra);
    _red  = MODIS_RED(can, spectra);
    _nir  = MODIS_NIR(can, spectra);

    return FT(2.5) * (_nir - _red) / (_nir + 6 * _red - FT(7.5) * _blue + 1)
);


"""

    MODIS_EVI2(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return EVI2 for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_EVI2 end;

MODIS_EVI2(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_EVI2(spac.CANOPY, config.SPECTRA);

MODIS_EVI2(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    _red = MODIS_RED(can, spectra);
    _nir = MODIS_NIR(can, spectra);

    return FT(2.5) * (_nir - _red) / (_nir + FT(2.4) * _red + 1)
);


"""

    MODIS_LSWI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return LSWI for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_LSWI end;

MODIS_LSWI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_LSWI(spac.CANOPY, config.SPECTRA);

MODIS_LSWI(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    _nir  = MODIS_NIR(can, spectra);
    _swir = read_spectrum(spectra.Λ, can.sensor_geometry.auxil.reflectance, FT(MODIS_BAND_7[1]), FT(MODIS_BAND_7[2]); steps=5);

    return (_nir - _swir) / (_nir + _swir)
);


"""

    MODIS_NDVI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return NDVI for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_NDVI end;

MODIS_NDVI(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_NDVI(spac.CANOPY, config.SPECTRA);

MODIS_NDVI(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    _red = MODIS_RED(can, spectra);
    _nir = MODIS_NIR(can, spectra);

    return (_nir - _red) / (_nir + _red)
);


"""

    MODIS_NIR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return near infrared band reflectance for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_NIR end;

MODIS_NIR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_NIR(spac.CANOPY, config.SPECTRA);

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

MODIS_NIRv(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_NIRv(spac.CANOPY, config.SPECTRA);

MODIS_NIRv(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    _red = MODIS_RED(can, spectra);
    _nir = MODIS_NIR(can, spectra);

    return (_nir - _red) / (_nir + _red) * _nir
);


"""

    MODIS_NIRvR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return NIRv radiance for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_NIRvR end;

MODIS_NIRvR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_NIRvR(spac.CANOPY, config.SPECTRA);

MODIS_NIRvR(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    _nir_rad = read_spectrum(spectra.Λ, can.sensor_geometry.auxil.e_sensor, FT(MODIS_BAND_2[1]), FT(MODIS_BAND_2[2]); steps=6);

    return MODIS_NDVI(can, spectra) * _nir_rad
);


"""

    MODIS_RED(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return red band reflectance for MODIS setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function MODIS_RED end;

MODIS_RED(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = MODIS_RED(spac.CANOPY, config.SPECTRA);

MODIS_RED(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ, can.sensor_geometry.auxil.reflectance, FT(MODIS_BAND_1[1]), FT(MODIS_BAND_1[2]); steps=6)
);


"""

    OCO2_SIF759(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 759 nm for OCO2 setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function OCO2_SIF759 end;

OCO2_SIF759(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = OCO2_SIF759(spac.CANOPY, config.SPECTRA);

OCO2_SIF759(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(OCO2_SIF_759[1]), FT(OCO2_SIF_759[2]); steps=4)
);


"""

    OCO2_SIF770(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 770 nm for OCO2 setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function OCO2_SIF770 end;

OCO2_SIF770(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = OCO2_SIF770(spac.CANOPY, config.SPECTRA);

OCO2_SIF770(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(OCO2_SIF_770[1]), FT(OCO2_SIF_770[2]); steps=4)
);


"""

    OCO3_SIF759(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 759 nm for OCO3 setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function OCO3_SIF759 end;

OCO3_SIF759(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = OCO3_SIF759(spac.CANOPY, config.SPECTRA);

OCO3_SIF759(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(OCO3_SIF_759[1]), FT(OCO3_SIF_759[2]); steps=4)
);


"""

    OCO3_SIF770(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 770 nm for OCO3 setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function OCO3_SIF770 end;

OCO3_SIF770(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = OCO3_SIF770(spac.CANOPY, config.SPECTRA);

OCO3_SIF770(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(OCO3_SIF_770[1]), FT(OCO3_SIF_770[2]); steps=4)
);


"""

    TROPOMI_SIF683(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 682.5 nm for TROPOMI setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function TROPOMI_SIF683 end;

TROPOMI_SIF683(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF683(spac.CANOPY, config.SPECTRA);

TROPOMI_SIF683(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(TROPOMI_SIF_683[1]), FT(TROPOMI_SIF_683[2]); steps=5)
);


"""

    TROPOMI_SIF740(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 740 nm for TROPOMI setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function TROPOMI_SIF740 end;

TROPOMI_SIF740(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF740(spac.CANOPY, config.SPECTRA);

TROPOMI_SIF740(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(740))
);


"""

    TROPOMI_SIF747(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 746.5 nm for TROPOMI setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function TROPOMI_SIF747 end;

TROPOMI_SIF747(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF747(spac.CANOPY, config.SPECTRA);

TROPOMI_SIF747(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(TROPOMI_SIF_747[1]), FT(TROPOMI_SIF_747[2]); steps=8)
);


"""

    TROPOMI_SIF751(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 750.5 nm for TROPOMI setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function TROPOMI_SIF751 end;

TROPOMI_SIF751(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF751(spac.CANOPY, config.SPECTRA);

TROPOMI_SIF751(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(TROPOMI_SIF_751[1]), FT(TROPOMI_SIF_751[2]); steps=5)
);
