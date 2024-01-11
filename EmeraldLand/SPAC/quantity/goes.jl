const GOES_R_ABI_BAND1 = [450, 490];    # BLUE
const GOES_R_ABI_BAND2 = [590, 690];    # RED
const GOES_R_ABI_BAND3 = [850, 890];    # NIR
const GOES_R_ABI_BAND4 = [1370, 1390];  # SWIR
const GOES_R_ABI_BAND5 = [1580, 1640];  # SWIR
const GOES_R_ABI_BAND6 = [2230, 2280];  # SWIR
const GOES_R_ABI_BANDS = [GOES_R_ABI_BAND1, GOES_R_ABI_BAND2, GOES_R_ABI_BAND3, GOES_R_ABI_BAND4, GOES_R_ABI_BAND5, GOES_R_ABI_BAND6];


#######################################################################################################################################################################################################
#
# Changes to these functions
# General
#     2024-Jan-11: add functions for GOES-R ABI geostaionary satellite
#
#######################################################################################################################################################################################################
"""

    GOES_R_ABI_BRFX(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, x::Int) where {FT}

Return band x reflectance for GOES-R ABI geostaionary satellite, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` type SPAC
- `x` band number

"""
function GOES_R_BRFX end;

GOES_R_BRFX(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, x::Int) where {FT} = GOES_R_BRFX(spac.canopy, config.SPECTRA, x);

GOES_R_BRFX(can::MultiLayerCanopy{FT}, spectra::ReferenceSpectra{FT}, x::Int) where {FT} = (
    return read_spectrum(spectra.Λ, can.sensor_geometry.auxil.reflectance, FT(GOES_R_ABI_BANDS[x][1]), FT(GOES_R_ABI_BANDS[x][2]); steps=4);
);

GOES_R_BRF1(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = GOES_R_BRFX(config, spac, 1);

GOES_R_BRF2(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = GOES_R_BRFX(config, spac, 2);

GOES_R_BRF3(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = GOES_R_BRFX(config, spac, 3);

GOES_R_BRF4(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = GOES_R_BRFX(config, spac, 4);

GOES_R_BRF5(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = GOES_R_BRFX(config, spac, 5);

GOES_R_BRF6(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = GOES_R_BRFX(config, spac, 6);
