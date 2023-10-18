const OCO2_SIF_759 = [758.17, 759.20];  # SIF 757
const OCO2_SIF_770 = [769.62, 770.28];  # SIF 770
const OCO3_SIF_759 = [758.26, 759.28];  # SIF 757
const OCO3_SIF_770 = [769.67, 770.34];  # SIF 770


#######################################################################################################################################################################################################
#
# Changes to these functions
# General
#     2022-Jun-13: add function to compute OCO2 SIF @ 758.7 and 770.0 nm
#     2022-Jun-13: add function to compute OCO3 SIF @ 758.8 and 770.0 nm
#
#######################################################################################################################################################################################################
"""

    OCO2_SIF759(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 759 nm for OCO2 setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function OCO2_SIF759 end;

OCO2_SIF759(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = OCO2_SIF759(spac.canopy, config.SPECTRA);

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

OCO2_SIF770(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = OCO2_SIF770(spac.canopy, config.SPECTRA);

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

OCO3_SIF759(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = OCO3_SIF759(spac.canopy, config.SPECTRA);

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

OCO3_SIF770(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = OCO3_SIF770(spac.canopy, config.SPECTRA);

OCO3_SIF770(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(OCO3_SIF_770[1]), FT(OCO3_SIF_770[2]); steps=4)
);
