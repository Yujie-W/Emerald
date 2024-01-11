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
    end;
end;
