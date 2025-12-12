"""

    site_driver_tuple(gmd::Union{Dict,OrderedDict}, wd::Union{Dict{String,Vector{FT}},OrderedDict{String,Vector{FT}}}) where {FT}

Prepare the Tuple of weather and trait drivers to drive the simulations, given
- `gmd` Dictionary of GriddingMachine data in a grid
- `wd` Dictionary of weather driver data in a grid

"""
function site_driver_tuple(gmd::Union{Dict,OrderedDict}, wd::Union{Dict{String,Vector{FT}},OrderedDict{String,Vector{FT}}}) where {FT}
    wd["B6F"    ] = resample(FT.(gmd["B6F"        ]), "1H", gmd["YEAR"]);
    wd["CO2"    ] = resample(FT.(gmd["CO2"        ]), "1H", gmd["YEAR"]);
    wd["CHL"    ] = resample(FT.(gmd["CHLOROPHYLL"]), "1H", gmd["YEAR"]);
    wd["CI"     ] = resample(FT.(gmd["CLUMPING"   ]), "1H", gmd["YEAR"]);
    wd["JMAX25" ] = resample(FT.(gmd["JMAX25"     ]), "1H", gmd["YEAR"]);
    wd["LAI"    ] = resample(FT.(gmd["LAI"        ]), "1H", gmd["YEAR"]);
    wd["VCMAX25"] = resample(FT.(gmd["VCMAX25"    ]), "1H", gmd["YEAR"]);

    # convert the DataFrame to NamedTuple
    return NamedTuple{Tuple(Symbol.(keys(wd)))}(Tuple([wd[k] for k in keys(wd)]))
end;
