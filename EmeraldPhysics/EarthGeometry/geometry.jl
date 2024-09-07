#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2024-Jan-11: add function to calculate the position of an object on the Earth's surface
#
#######################################################################################################################################################################################################
"""

    object_position(lat::FT, lon::FT, h::FT) where {FT}

Return the position of an object on the Earth's surface (x for center to [0,0], y for center to [0,90], z for center to north pole), given
- `lat` Latitude in °
- `lon` Longitude in °
- `h` Height in m

"""
function object_position(lat::FT, lon::FT, h::FT) where {FT}
    r_eq = R_EQUATOR(FT);
    r_pl = R_POLAR(FT);

    f = 1 - r_pl / r_eq;
    c = 1 / sqrt(1 + f * (f - 2) * sind(lat)^2);
    sq = c * (1 - f) ^ 2;

    achcp = (r_eq * c + h) * cosd(lat);
    x = achcp * cosd(lon);
    y = achcp * sind(lon);
    z = (r_eq * sq + h) * sind(lat);

    return x, y, z
end;


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2024-Jan-11: add function to calculate the viewer angles of an object on the Earth's surface
# Bug fixes
#     2024-Sep-06: fix the calculation of the viewer zenith angle (was using elevation angle)
#
#######################################################################################################################################################################################################
"""

    viewer_angles(sat_lat::FT, sat_lon::FT, sat_h::FT, lat::FT, lon::FT, h::FT = FT(0)) where {FT}

Return the viewer zenith and azimuth angles, given
- `sat_lat` Satellite latitude in °
- `sat_lon` Satellite longitude in °
- `sat_h` Satellite height in m
- `lat` Target latitude in °
- `lon` Target longitude in °
- `h` Target height on Earth surface in m

"""
function viewer_angles(sat_lat::FT, sat_lon::FT, sat_h::FT, lat::FT, lon::FT, h::FT = FT(0)) where {FT}
    (sat_x, sat_y, sat_z) = object_position(sat_lat, sat_lon, sat_h);
    (tar_x, tar_y, tar_z) = object_position(lat, lon, h);

    rx = sat_x - tar_x;
    ry = sat_y - tar_y;
    rz = sat_z - tar_z;

    top_s = sind(lat) * cosd(lon) * rx + sind(lat) * sind(lon) * ry - cosd(lat) * rz;
    top_e = -sind(lon) * rx + cosd(lon) * ry
    top_z = cosd(lat) * cosd(lon) * rx + cosd(lat) * sind(lon) * ry + sind(lat) * rz;

    # viewer azimuth angle
    vaa = atand(-top_e, top_s) + 180;

    # viewer zenith angle
    r_g = sqrt(rx^2 + ry^2 + rz^2);
    vza = 90 - asind(min(1, top_z / r_g));

    return vza, vaa
end;
