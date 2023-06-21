#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-May-04: add new function animate_data!
#     2023-May-17: add new method for 2D data
#
#######################################################################################################################################################################################################
"""

    animate_data!(
                xs::Vector{<:Number},
                ys::Vector{<:Number},
                data::Array{<:Number,3};
                cmap::Symbol = :viridis,
                filename::Union{Nothing,String} = nothing,
                figsize::Tuple{Int,Int} = (400,300),
                fps::Int = 15,
                ismap::Bool = true,
                titles::Vector{String} = ["" for _ in axes(data,3)],
                vmin_vmax::Union{Nothing,Tuple{Number,Number}} = (minimum(data), maximum(data)),
                xorg::Bool = false)
    animate_data!(
                xs::Vector{<:Number},
                ys::Vector{<:Number},
                data::Matrix{<:Number};
                cmap::Symbol = :viridis,
                filename::Union{Nothing,String} = nothing,
                figsize::Tuple{Int,Int} = (400,300),
                ismap::Bool = true,
                vmin_vmax::Union{Nothing,Tuple{Number,Number}} = (nanmin(data), nanmax(data)),
                xorg::Bool = false)

Animate the netcdf file, given
- `xs` X-axis (first dimension of data)
- `ys` Y-axis (second dimension of data)
- `data` 3D or 2D data
- `cmap` Color scheme to use
- `filename` File name or path to save the animation (to gif or mp4; default is nothing)
- `figsize` Cavas size in pixels
- `fps` Frame per second
- `ismap` Whether the data is a map
- `titles` Title of each frame
- `vmin_vmax` Min and max value of the colorbar
- `xorg` If false, set the `ENV["GKSwstype"]` to 100 to diable warnings

"""
function animate_data! end

animate_data!(
            xs::Vector{<:Number},
            ys::Vector{<:Number},
            data::Array{<:Number,3};
            cmap::Symbol = :viridis,
            filename::Union{Nothing,String} = nothing,
            figsize::Tuple{Int,Int} = (400,300),
            fps::Int = 15,
            ismap::Bool = true,
            titles::Vector{String} = ["" for _ in axes(data,3)],
            vmin_vmax::Union{Nothing,Tuple{Number,Number}} = (nanmin(data), nanmax(data)),
            xorg::Bool = false) = (
    # disable GKS warning if xorg is false
    xorg ? nothing : ENV["GKSwstype"] = "100";

    _make_frame(i) = (
        _fig = heatmap(xs, ys, data[:,:,i]'; c = cmap, clim = vmin_vmax, size = figsize, title = titles[i]);
        ismap ? set_style!(_fig, GlobalMapStyle(); latitude_360 = (maximum(xs) > 180)) : nothing;
    );

    return animation(_make_frame, collect(axes(data,3)); filename = filename, fps = fps)
);

animate_data!(
            xs::Vector{<:Number},
            ys::Vector{<:Number},
            data::Matrix{<:Number};
            cmap::Symbol = :viridis,
            filename::Union{Nothing,String} = nothing,
            figsize::Tuple{Int,Int} = (400,300),
            ismap::Bool = true,
            vmin_vmax::Union{Nothing,Tuple{Number,Number}} = (nanmin(data), nanmax(data)),
            xorg::Bool = false) = (
    # disable GKS warning if xorg is false
    xorg ? nothing : ENV["GKSwstype"] = "100";

    _fig = heatmap(xs, ys, data[:,:]'; c = cmap, clim = vmin_vmax, size = figsize);
    ismap ? set_style!(_fig, GlobalMapStyle(); latitude_360 = (maximum(xs) > 180)) : nothing;

    return isnothing(filename) ? _fig : savefig(_fig, filename)
);
