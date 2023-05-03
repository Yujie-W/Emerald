#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Feb-15: add new function
#     2023-May-03: add new method to plot data directly
#
#######################################################################################################################################################################################################
"""

    animate_nc!(datafile::String, varname::String; filename::Union{Nothing,String} = nothing, fps::Int = 15)
    animate_nc!(lats::Vector{<:Number}, lons::Vector{<:Number}, data::Array{<:Number,3}; filename::Union{Nothing,String} = nothing, fps::Int = 15)

Animate the netcdf file, given
- `datafile` Path to netcdf dataset
- `varname` Variable name
- `filename` File name or path to save the animation (to gif or mp4; default is nothing)
- `fps` Frame per second (default is 15)
- `lats` Latitude (y axis)
- `lons` Longitude (x axis)
- `data` 3D data (x, y, z)

"""
function animate_nc! end

animate_nc!(datafile::String, varname::String; filename::Union{Nothing,String} = nothing, fps::Int = 15) = (
    _lats = read_nc(datafile, "lat");
    _lons = read_nc(datafile, "lon");
    _data = read_nc(datafile, varname);

    return animate_nc!(_lats, _lons, _data; filename = filename, fps = fps)
);

animate_nc!(lats::Vector{<:Number}, lons::Vector{<:Number}, data::Array{<:Number,3}; filename::Union{Nothing,String} = nothing, fps::Int = 15) = (
    _make_frame(i) = (
        _tmp = data[:,:,i]';
        _fig = heatmap(lons, lats, _tmp);
        set_style!(_fig, GlobalMapStyle(); latitude_360 = (maximum(lons) > 180));
    );

    return animation(_make_frame, collect(axes(data,3)); filename = filename, fps = fps)
);
