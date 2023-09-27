module EarthGeometry

using ..Constant: YEAR_D


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Apr-26: add function daytime_length()
#
#######################################################################################################################################################################################################
"""

    day_length(lat::FT, fdoy::FT) where {FT}

Return the daytime length, given
- `lat` Latitude in °
- `fdoy` Day of year (digits after decimal for time of day)

"""
function day_length(lat::FT, fdoy::FT) where {FT}
    _deg  = 360 / YEAR_D(FT) * (fdoy + 10);
    _decl = -FT(23.44) * cosd(_deg);

    if -tand(lat) * tand(_decl) <= -1
        return FT(24)
    elseif -tand(lat) * tand(_decl) > 1
        return FT(0)
    else
        return 2 * acosd(-tand(lat) * tand(_decl)) / 15
    end
end


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Jun-09: add function solar_azimuth_angle
#
#######################################################################################################################################################################################################
"""

    solar_azimuth_angle(lat::FT, decl::FT, lha::FT) where {FT}
    solar_azimuth_angle(lat::FT, day::FT, hour::FT, minute::FT) where {FT}
    solar_azimuth_angle(lat::FT, fdoy::FT) where {FT}

Return the solar azimuth angle, given
- `lat` Latitude in °
- `decl` Declination of the Sun in °
- `lha` Local hour angle in °
- `day` Day of year
- `hour` Hour of day
- `minute` Minute of hour
- `fdoy` Day of year (digits after decimal for time of day)

"""
function solar_azimuth_angle end

solar_azimuth_angle(lat::FT, decl::FT, lha::FT) where {FT} = (
    _cos_sza = sind(lat) * sind(decl) + cosd(lat) * cosd(decl) * cosd(lha);

    if _cos_sza == 1
        return lat >= 0 ? FT(180) : FT(0)
    end;

    _sin_sza = sqrt(1 - _cos_sza^2);
    _cos_saa = (sind(decl) - sind(lat) * _cos_sza) / (cosd(lat) * _sin_sza);

    return lha >= 0 ? 360 - acosd(_cos_saa) : acosd(_cos_saa)
);

solar_azimuth_angle(lat::FT, day::FT, hour::FT, minute::FT) where {FT} = (
    _decl = solar_declination_angle(day + (hour + minute / 60) / 24);
    _lha  = (hour - 12) * 15;

    return solar_azimuth_angle(lat, _decl, _lha)
);

solar_azimuth_angle(lat::FT, fdoy::FT) where {FT} = (
    _decl = solar_declination_angle(fdoy);
    _lha  = ((fdoy % 1) - FT(0.5)) * 360;

    return solar_azimuth_angle(lat, _decl, _lha)
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Apr-26: add function solar_declination_angle
#
#######################################################################################################################################################################################################
"""

    solar_declination_angle(fdoy::FT) where {FT}

Return the solar declination angle, given
- `fdoy` Day of year

"""
function solar_declination_angle(fdoy::FT) where {FT}
    return -FT(23.44) * cosd(360 / YEAR_D(FT) * (fdoy + 10))
end


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Sep-08: move function from SoilPlantAirContinuum.jl
#     2022-Sep-08: make YEAR_D an option so that one can use ClimaCache value to replace it
#     2022-Oct-19: use YEAR_D function from EmeraldConstants
#
#######################################################################################################################################################################################################
"""

    solar_zenith_angle(lat::FT, decl::FT, lha::FT) where {FT}
    solar_zenith_angle(lat::FT, day::FT, hour::FT, minute::FT) where {FT}
    solar_zenith_angle(lat::FT, fdoy::FT) where {FT}

Return the solar zenith angle, given
- `lat` Latitude in °
- `decl` Declination of the Sun in °
- `lha` Local hour angle in °
- `day` Day of year
- `hour` Hour of day
- `minute` Minute of hour
- `fdoy` Day of year (digits after decimal for time of day)

"""
function solar_zenith_angle end

solar_zenith_angle(lat::FT, decl::FT, lha::FT) where {FT} = (
    _cos_sza = sind(lat) * sind(decl) + cosd(lat) * cosd(decl) * cosd(lha);

    return acosd(_cos_sza)
);

solar_zenith_angle(lat::FT, day::FT, hour::FT, minute::FT) where {FT} = (
    _decl = solar_declination_angle(day + (hour + minute / 60) / 24);
    _lha  = (hour - 12) * 15;

    return solar_zenith_angle(lat, _decl, _lha)
);

solar_zenith_angle(lat::FT, fdoy::FT) where {FT} = (
    _decl = solar_declination_angle(fdoy);
    _lha  = ((fdoy % 1) - FT(0.5)) * 360;

    return solar_zenith_angle(lat, _decl, _lha)
);


end # module
