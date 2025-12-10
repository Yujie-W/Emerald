const OCO2_SIF_759 = [758.17, 759.20];  # SIF 757
const OCO2_SIF_770 = [769.62, 770.28];  # SIF 770
const OCO3_SIF_759 = [758.26, 759.28];  # SIF 757
const OCO3_SIF_770 = [769.67, 770.34];  # SIF 770
const TROPOMI_SIF_683 = [680, 685];     # SIF 683
const TROPOMI_SIF_747 = [735, 758];     # SIF 747
const TROPOMI_SIF_751 = [743, 758];     # SIF 751


"""

    OCO2_SIF759(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 759 nm for OCO2 setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function OCO2_SIF759 end;

OCO2_SIF759(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = OCO2_SIF759(spac.canopy, config.CONSTANTS.SPECTRA);

OCO2_SIF759(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(OCO2_SIF_759[1]), FT(OCO2_SIF_759[2]); steps=4)
);


"""

    OCO2_SIF770(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 770 nm for OCO2 setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function OCO2_SIF770 end;

OCO2_SIF770(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = OCO2_SIF770(spac.canopy, config.CONSTANTS.SPECTRA);

OCO2_SIF770(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(OCO2_SIF_770[1]), FT(OCO2_SIF_770[2]); steps=4)
);


"""

    OCO3_SIF759(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 759 nm for OCO3 setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function OCO3_SIF759 end;

OCO3_SIF759(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = OCO3_SIF759(spac.canopy, config.CONSTANTS.SPECTRA);

OCO3_SIF759(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(OCO3_SIF_759[1]), FT(OCO3_SIF_759[2]); steps=4)
);


"""

    OCO3_SIF770(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 770 nm for OCO3 setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function OCO3_SIF770 end;

OCO3_SIF770(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = OCO3_SIF770(spac.canopy, config.CONSTANTS.SPECTRA);

OCO3_SIF770(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(OCO3_SIF_770[1]), FT(OCO3_SIF_770[2]); steps=4)
);


"""

    TROPOMI_SIF683(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 682.5 nm for TROPOMI setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function TROPOMI_SIF683 end;

TROPOMI_SIF683(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF683(spac.canopy, config.CONSTANTS.SPECTRA);

TROPOMI_SIF683(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(TROPOMI_SIF_683[1]), FT(TROPOMI_SIF_683[2]); steps=5)
);


"""

    TROPOMI_SIF740(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 740 nm for TROPOMI setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function TROPOMI_SIF740 end;

TROPOMI_SIF740(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF740(spac.canopy, config.CONSTANTS.SPECTRA);

TROPOMI_SIF740(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return interpolate_data(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(740))
);


"""

    TROPOMI_SIF747(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 746.5 nm for TROPOMI setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function TROPOMI_SIF747 end;

TROPOMI_SIF747(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF747(spac.canopy, config.CONSTANTS.SPECTRA);

TROPOMI_SIF747(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(TROPOMI_SIF_747[1]), FT(TROPOMI_SIF_747[2]); steps=8)
);


"""

    TROPOMI_SIF751(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return SIF @ 750.5 nm for TROPOMI setup, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC

"""
function TROPOMI_SIF751 end;

TROPOMI_SIF751(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = TROPOMI_SIF751(spac.canopy, config.CONSTANTS.SPECTRA);

TROPOMI_SIF751(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}) where {FT} = (
    return read_spectrum(spectra.Λ_SIF, can.sensor_geometry.auxil.sif_obs, FT(TROPOMI_SIF_751[1]), FT(TROPOMI_SIF_751[2]); steps=5)
);


"""

    ΣSIF(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return the total SIF at top of the canopy after reabsorption in W m⁻² per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function ΣSIF end;

ΣSIF(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = (
    (; SPECTRA) = config.CONSTANTS;
    sun_geo = spac.canopy.sun_geometry;

    return sun_geo.auxil.e_sifꜛ[:,1]' * SPECTRA.ΔΛ_SIF / 1000
);


"""

    ΣSIF_CHL(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return the total SIF at chloroplast level (without any reabsorption) in W m⁻² per ground area, given
- `config` `SPACConfig` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function ΣSIF_CHL end;

ΣSIF_CHL(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = (
    (; SPECTRA) = config.CONSTANTS;
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute SIF in energy unit before reabsorption within leaves (W m⁻²)
    Σsif::FT = 0;
    for i in eachindex(leaves)
        Σsif += view(canopy.sun_geometry.auxil.e_sif_chl,:,i)' * SPECTRA.ΔΛ_SIF / 1000;
    end;

    return Σsif
);


"""

    ΣSIF_LEAF(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return the total SIF at leaf level after reabsorption in W m⁻² per ground area, given
- `config` `SPACConfig` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function ΣSIF_LEAF end;

ΣSIF_LEAF(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = (
    (; SPECTRA) = config.CONSTANTS;
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute SIF in energy unit after reabsorption within leaves (W m⁻²)
    Σsif::FT = 0;
    for i in eachindex(leaves)
        Σsif += view(canopy.sun_geometry.auxil.e_sifꜜ_layer,:,i)' * SPECTRA.ΔΛ_SIF / 1000;
        Σsif += view(canopy.sun_geometry.auxil.e_sifꜛ_layer,:,i)' * SPECTRA.ΔΛ_SIF / 1000;
    end;

    return Σsif
);
