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

    day_length(lat::FT, fdoy::FT) where {FT<:AbstractFloat}

Return the daytime length, given
- `lat` Latitude in °
- `fdoy` Day of year (digits after decimal for time of day)

"""
function day_length(lat::FT, fdoy::FT) where {FT<:AbstractFloat}
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
#     2023-Apr-26: add function solar_declination_angle
#
#######################################################################################################################################################################################################
"""

    solar_declination_angle(fdoy::FT) where {FT<:AbstractFloat}

Return the solar declination angle, given
- `fdoy` Day of year

"""
function solar_declination_angle(fdoy::FT) where {FT<:AbstractFloat}
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

    solar_zenith_angle(lat::FT, decl::FT, lha::FT) where {FT<:AbstractFloat}
    solar_zenith_angle(lat::FT, day::FT, hour::FT, minute::FT) where {FT<:AbstractFloat}
    solar_zenith_angle(lat::FT, fdoy::FT) where {FT<:AbstractFloat}

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

solar_zenith_angle(lat::FT, decl::FT, lha::FT) where {FT<:AbstractFloat} = (
    _cos_sza = sind(lat) * sind(decl) + cosd(lat) * cosd(decl) * cosd(lha);

    return acosd(_cos_sza)
);

solar_zenith_angle(lat::FT, day::FT, hour::FT, minute::FT) where {FT<:AbstractFloat} = (
    _decl = solar_declination_angle(day + (hour + minute / 60) / 24);
    _lha  = (hour - 12) * 15;

    return solar_zenith_angle(lat, _decl, _lha)
);

solar_zenith_angle(lat::FT, fdoy::FT) where {FT<:AbstractFloat} = (
    _decl = solar_declination_angle(fdoy);
    _lha  = ((fdoy % 1) - FT(0.5)) * 360;

    return solar_zenith_angle(lat, _decl, _lha)
);


end # module
