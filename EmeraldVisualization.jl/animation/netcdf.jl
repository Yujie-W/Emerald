#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Feb-15: add new function
#
#######################################################################################################################################################################################################
"""

    animate_nc!(datafile::String, varname::String; kwargs...)

Animate the netcdf file, given
- `datafile` Path to netcdf dataset
- `varname` Variable name
- `kwargs` Keyword arguments of [`animate_data!`](@ref)

"""
function animate_nc!(datafile::String, varname::String; kwargs...)
    _lats = read_nc(datafile, "lat");
    _lons = read_nc(datafile, "lon");
    _data = read_nc(datafile, varname);

    return animate_data!(_lons, _lats, _data; kwargs...)
end
