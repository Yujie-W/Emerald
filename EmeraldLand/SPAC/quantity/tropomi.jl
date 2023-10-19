const TROPOMI_SIF_683 = [680, 685];     # SIF 683
const TROPOMI_SIF_747 = [735, 758];     # SIF 747
const TROPOMI_SIF_751 = [743, 758];     # SIF 751


#######################################################################################################################################################################################################
#
# Changes to these functions
# General
#     2022-Jun-13: add function to compute TROPOMI SIF @ 682.5, 740.0, 746.5, and 750.5 nm
#
#######################################################################################################################################################################################################
"""

    TROPOMI_SIF683(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 682.5 nm for TROPOMI setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function TROPOMI_SIF683 end;

TROPOMI_SIF683(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF683(spac.canopy, config.SPECTRA);

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

TROPOMI_SIF740(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF740(spac.canopy, config.SPECTRA);

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

TROPOMI_SIF747(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF747(spac.canopy, config.SPECTRA);

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

TROPOMI_SIF751(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF751(spac.canopy, config.SPECTRA);

TROPOMI_SIF751(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(TROPOMI_SIF_751[1]), FT(TROPOMI_SIF_751[2]); steps=5)
);
