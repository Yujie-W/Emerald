#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jan-20: refactor the function
#
#######################################################################################################################################################################################################
"""

    plot_ellipse!(
                ax::PyObject,
                xy::Tuple{Number, Number};
                alpha::Number = 0.5,
                angle::Number = 0,
                color::String = "black",
                edgecolor::String = color,
                facecolor::String = color,
                height::Number = 10,
                width::Number = 10)

Plot an ellipse on axis, given
- `ax` Axis to plot on
- `xy` Center of the ellipse
- `alpha` Transparency of the ellipse
- `angle` Rotation angle of the ellipse
- `color` Color of the ellipse
- `edgecolor` Edgecolor of the ellipse
- `facecolor` Face color of the ellipse
- `height` Height of the ellipse
- `width` Width of the ellipse

"""
function plot_ellipse!(
            ax::PyObject,
            xy::Tuple{Number, Number};
            alpha::Number = 0.5,
            angle::Number = 0,
            color::String = "black",
            edgecolor::String = color,
            facecolor::String = color,
            height::Number = 10,
            width::Number = 10)
    _outl = PATCHES.Ellipse(; xy = xy, width = width, height = height, edgecolor = edgecolor, facecolor = facecolor, angle = angle, clip_on = false, alpha = alpha);
    ax.add_artist(_outl);

    return nothing
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jan-20: refactor the function
#
#######################################################################################################################################################################################################
"""

    plot_stoma!(ax::PyObject, xy::Tuple{Number, Number}; angle::Number = 0, height::Number = 10, stoma::Number = 0.2, width::Number = 10)

Plot a stoma on the axis, given
- `ax` Axis to plot on
- `xy` Center of the stoma
- `angle` Rotation angle of the stoma
- `height` Height of the stoma
- `stoma` Stomatal pore width ratio
- `width` Width of the stoma

"""
function plot_stoma!(ax::PyObject, xy::Tuple{Number, Number}; angle::Number = 0, height::Number = 10, stoma::Number = 0.2, width::Number = 10)
    plot_ellipse!(ax, xy; angle = angle, color = "green", height = height, width = width);
    plot_ellipse!(ax, xy; angle = angle, color = "white", height = height * 0.9, width = width * stoma);

    return nothing
end
