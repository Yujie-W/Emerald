#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Dec-07: refactor the function to canvas
#     2023-May-17: add option polar
#
#######################################################################################################################################################################################################
"""

    canvas(id::Union{Int, String};
           axids::Union{Nothing, Vector{Int}} = nothing,
           dpi::Int = 100,
           ncol::Int = 1,
           nrow::Int = 1,
           polar::Union{Bool, Vector{Bool}} = false,
           figsize::Tuple = (0.5 + ncol * 3, 0.5 + nrow * 3)

Create the canvas, given
- `id` ID of the figure
- `axids` Given indicies of the subplots in the figure
- `dpi` Given pixels per inch
- `ncol` Number of columns in the figure
- `nrow` Number of rows in the figure
- `polar` Whether the panels are polar projected
- `figsize` Given canvas size

"""
function canvas(
            id::Union{Int, String};
            axids::Union{Nothing, Vector{Int}} = nothing,
            dpi::Int = 100,
            ncol::Int = 1,
            nrow::Int = 1,
            polar::Union{Bool, Vector{Bool}} = false,
            figsize::Tuple = (0.5 + ncol * 3, 0.5 + nrow * 3))
    # create a clean figure
    _fig = figure(id, figsize = figsize, dpi = dpi);
    _fig.clear();

    # if axids is empty
    _axids = isnothing(axids) ? collect(1:nrow*ncol) : axids;

    # add axes
    for _i in eachindex(_axids)
        if polar isa Bool
            _fig.add_subplot(nrow, ncol, _axids[_i], polar = polar);
        else
            _fig.add_subplot(nrow, ncol, _axids[_i], polar = polar[_i]);
        end;
    end;

    return _fig
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Dec-07: refactor the function to canvas
#
#######################################################################################################################################################################################################
"""

    save_canvas!(fig::Figure, figname::String; folder::String = "", formats::Vector{String} = ["pdf"], tight::Bool = true)

Save the figure, given
- `fig` PyPlot Figure (canvas)
- `figname` Figure name to save
- `folder` Containing folder
- `formats` Vector of formats, must be within [jpg, pdf, png]
- `tight` If true, use tight layout

"""
function save_canvas!(fig::Figure, figname::String; folder::String = "./", formats::Vector{String} = ["pdf"], tight::Bool = true)
    if tight fig.set_tight_layout(true); end;

    for _format in formats
        @assert _format in ["jpg", "pdf", "png"] "Output format must be one of jpg, pdf, and png!";
        fig.savefig("$(folder)/$(figname).$(_format)", bbox_inches="tight");
    end;

    return nothing
end
