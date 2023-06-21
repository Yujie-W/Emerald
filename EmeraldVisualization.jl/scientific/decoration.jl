#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Dec-07: refactor the functions that set title and axes to decorate!
#     2022-Dec-07: add options to set titles
#     2022-Dec-07: add options to set x and y limits,
#     2022-Dec-07: add options to set x and y ticks and ticklabels
#     2022-Dec-07: add options to set x and y labels
#
#######################################################################################################################################################################################################
"""

`decorate!` supports two methods: one to decorate a vector of axes and one to decorate a single axis.

"""
function decorate! end

"""

    decorate!(axs::Vector;
              add_title::Bool = true,
              title_capital::Bool = false,
              title_dot::Bool = false,
              title_fontsize::Int = 16,
              title_labels::Union{Nothing, Vector{String}} = nothing,
              title_loc::String = "left",
              title_number::Bool = false,
              title_parentheses::Bool = true,
              use_latex::Bool = false,
              xaxis_label_fontsize::Int = 16,
              xaxis_labels::Union{Nothing, String, Vector{String}} = nothing,
              xaxis_lims::Union{Nothing, Tuple{Number,Number}, Vector} = nothing,
              xaxis_ticklabels::Union{Nothing, Vector} = nothing,
              xaxis_ticks::Union{Nothing, Vector} = nothing,
              yaxis_label_fontsize::Int = 16,
              yaxis_labels::Union{Nothing, String, Vector{String}} = nothing,
              yaxis_lims::Union{Nothing, Tuple{Number,Number}, Vector} = nothing,
              yaxis_ticklabels::Union{Nothing, Vector} = nothing,
              yaxis_ticks::Union{Nothing, Vector} = nothing)

Decorate the axes, given
- `axes` Vector of PyPlot axis
- `add_title` If true, add title to each axis
- `title_capital` If true, use capital letter in axis title
- `title_dot` If true, add a dot between title and extra label
- `title_fontsize` Title font size
- `title_labels` Vector of extra labels
- `title_loc` Title location
- `title_number` If true, use number instead of letter as the index
- `title_parentheses` If true, add `()` along with the index
- `use_latex` Use LaTeX to render the index
- `xaxis_label_fontsize` X-axis label font size
- `xaxis_labels` X-axis labels
- `xaxis_lims` X-axis limits
- `xaxis_ticklabels` X-axis tick labels
- `xaxis_ticks` X-axis ticks
- `yaxis_label_fontsize` Y-axis label font size
- `yaxis_labels` Y-axis labels
- `yaxis_lims` Y-axis limits
- `yaxis_ticklabels` Y-axis tick labels
- `yaxis_ticks` Y-axis ticks

"""
decorate!(axs::Vector;
          add_title::Bool = true,
          title_capital::Bool = false,
          title_dot::Bool = false,
          title_fontsize::Int = 16,
          title_labels::Union{Nothing, Vector{String}} = nothing,
          title_loc::String = "left",
          title_number::Bool = false,
          title_parentheses::Bool = true,
          use_latex::Bool = false,
          xaxis_label_fontsize::Int = 16,
          xaxis_labels::Union{Nothing, String, Vector{String}} = nothing,
          xaxis_lims::Union{Nothing, Tuple{Number,Number}, Vector} = nothing,
          xaxis_ticklabels::Union{Nothing, Vector} = nothing,
          xaxis_ticks::Union{Nothing, Vector} = nothing,
          yaxis_label_fontsize::Int = 16,
          yaxis_labels::Union{Nothing, String, Vector{String}} = nothing,
          yaxis_lims::Union{Nothing, Tuple{Number,Number}, Vector} = nothing,
          yaxis_ticklabels::Union{Nothing, Vector} = nothing,
          yaxis_ticks::Union{Nothing, Vector} = nothing) = (
    # make sure the dimenssions agree
    if !isnothing(title_labels)
        @assert length(axs) == length(title_labels) "Provided axes and title labels do not match in length!";
    end;
    if typeof(xaxis_labels) <: Vector
        @assert length(axs) == length(xaxis_labels) "Provided axes and x-axis labels do not match in length!";
    end;
    if typeof(yaxis_labels) <: Vector
        @assert length(axs) == length(yaxis_labels) "Provided axes and x-axis labels do not match in length!";
    end;
    if typeof(xaxis_lims) <: Vector
        @assert length(axs) == length(xaxis_lims) "Provided axes and x-axis limits do not match in length!";
    end;
    if typeof(yaxis_lims) <: Vector
        @assert length(axs) == length(yaxis_lims) "Provided axes and y-axis limits do not match in length!";
    end;
    if eltype(xaxis_ticks) <: Vector
        @assert length(axs) == length(xaxis_ticks) "Provided axes and x-axis ticks do not match in length!";
    end;
    if eltype(yaxis_ticks) <: Vector
        @assert length(axs) == length(yaxis_ticks) "Provided axes and y-axis ticks do not match in length!";
    end;
    if eltype(xaxis_ticklabels) <: Vector
        @assert length(axs) == length(xaxis_ticklabels) "Provided axes and x-axis ticklabels do not match in length!";
    end;
    if eltype(yaxis_ticklabels) <: Vector
        @assert length(axs) == length(yaxis_ticklabels) "Provided axes and y-axis ticklabels do not match in length!";
    end;

    # decorate title for each panel
    if add_title
        LETTERS = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"];
        NUMBERS = [string(_i) for _i in 1:100];
        _labels = title_number ? NUMBERS : LETTERS;
        for _id in eachindex(axs)
            _title = "";
            _label = _labels[_id];

            use_latex                ? _title *= "\\textbf{"             : nothing;
            title_parentheses        ? _title *= "("                     : nothing;
            title_capital            ? _title *= uppercase(_label)       : _title *= _label;
            title_parentheses        ? _title *= ")"                     : nothing;
            title_dot                ? _title *= "."                     : nothing;
            use_latex                ? _title *= "}"                     : nothing;
            !isnothing(title_labels) ? _title *= " " * title_labels[_id] : nothing;

            decorate!(axs[_id]; title = _title, title_fontsize = title_fontsize, title_loc = title_loc);
        end;
    end;

    # decorate x-axis and y-axis for each panel
    if !isnothing(xaxis_labels)
        for _id in eachindex(axs)
            _axis_label = typeof(xaxis_labels) <: Vector ? xaxis_labels[_id] : xaxis_labels;
            decorate!(axs[_id]; xaxis_label = _axis_label, xaxis_label_fontsize = xaxis_label_fontsize);
        end;
    end;
    if !isnothing(yaxis_labels)
        for _id in eachindex(axs)
            _axis_label = typeof(yaxis_labels) <: Vector ? yaxis_labels[_id] : yaxis_labels;
            decorate!(axs[_id]; yaxis_label = _axis_label, yaxis_label_fontsize = yaxis_label_fontsize);
        end;
    end;
    if !isnothing(xaxis_lims)
        for _id in eachindex(axs)
            typeof(xaxis_lims) <: Vector ? decorate!(axs[_id]; xaxis_lims = xaxis_lims[_id]) : decorate!(axs[_id]; xaxis_lims = xaxis_lims);
        end;
    end;
    if !isnothing(yaxis_lims)
        for _id in eachindex(axs)
            typeof(yaxis_lims) <: Vector ? decorate!(axs[_id]; yaxis_lims = yaxis_lims[_id]) : decorate!(axs[_id]; yaxis_lims = yaxis_lims);
        end;
    end;
    if !isnothing(xaxis_ticks)
        for _id in eachindex(axs)
            eltype(xaxis_ticks) <: Vector ? decorate!(axs[_id]; xaxis_ticks = xaxis_ticks[_id]) : decorate!(axs[_id]; xaxis_ticks = xaxis_ticks);
        end;
    end;
    if !isnothing(yaxis_ticks)
        for _id in eachindex(axs)
            eltype(yaxis_ticks) <: Vector ? decorate!(axs[_id]; yaxis_ticks = yaxis_ticks[_id]) : decorate!(axs[_id]; yaxis_ticks = yaxis_ticks);
        end;
    end;
    if !isnothing(xaxis_ticklabels)
        for _id in eachindex(axs)
            eltype(xaxis_ticklabels) <: Vector ? decorate!(axs[_id]; xaxis_ticklabels = xaxis_ticklabels[_id]) : decorate!(axs[_id]; xaxis_ticklabels = xaxis_ticklabels);
        end;
    end;
    if !isnothing(yaxis_ticklabels)
        for _id in eachindex(axs)
            eltype(yaxis_ticklabels) <: Vector ? decorate!(axs[_id]; yaxis_ticklabels = yaxis_ticklabels[_id]) : decorate!(axs[_id]; yaxis_ticklabels = yaxis_ticklabels);
        end;
    end;

    return nothing
);

"""

    decorate!(ax;
              title::Union{Nothing, String} = nothing,
              title_fontsize::Int = 16,
              title_loc::String = "left",
              use_latex::Bool = false,
              xaxis_label::Union{Nothing, String} = nothing,
              xaxis_label_fontsize::Int = 16,
              xaxis_lims::Union{Nothing, Tuple{Number,Number}} = nothing,
              xaxis_ticklabels::Union{Nothing, Vector} = nothing,
              xaxis_ticks::Union{Nothing, Vector} = nothing,
              yaxis_label::Union{Nothing, String} = nothing,
              yaxis_label_fontsize::Int = 16,
              yaxis_lims::Union{Nothing, Tuple{Number,Number}} = nothing,
              yaxis_ticklabels::Union{Nothing, Vector} = nothing,
              yaxis_ticks::Union{Nothing, Vector} = nothing)

Decorate a single axis, given
- `ax` PyPlot axis
- `title` Axis title
- `title_fontsize` Title font size
- `title_loc` Title location
- `use_latex` Use LaTeX to render the index
- `xaxis_label` X-axis label
- `xaxis_label_fontsize` X-axis label font size
- `xaxis_lims` X-axis limits
- `xaxis_ticklabels` X-axis tick labels
- `xaxis_ticks` X-axis ticks
- `yaxis_label` Y-axis label
- `yaxis_label_fontsize` Y-axis label font size
- `yaxis_lims` Y-axis limits
- `yaxis_ticklabels` Y-axis tick labels
- `yaxis_ticks` Y-axis ticks

"""
decorate!(ax;
          title::Union{Nothing, String} = nothing,
          title_fontsize::Int = 16,
          title_loc::String = "left",
          use_latex::Bool = false,
          xaxis_label::Union{Nothing, String} = nothing,
          xaxis_label_fontsize::Int = 16,
          xaxis_lims::Union{Nothing, Tuple{Number,Number}} = nothing,
          xaxis_ticklabels::Union{Nothing, Vector} = nothing,
          xaxis_ticks::Union{Nothing, Vector} = nothing,
          yaxis_label::Union{Nothing, String} = nothing,
          yaxis_label_fontsize::Int = 16,
          yaxis_lims::Union{Nothing, Tuple{Number,Number}} = nothing,
          yaxis_ticklabels::Union{Nothing, Vector} = nothing,
          yaxis_ticks::Union{Nothing, Vector} = nothing) = (
    # decorate title
    if !isnothing(title)
        if use_latex
            ax.set_title(title, fontsize = title_fontsize, loc = title_loc);
        else
            ax.set_title(title, fontsize = title_fontsize, loc = title_loc, fontweight="bold");
        end;
    end;

    # decorate x-axis and y-axis
    if !isnothing(xaxis_label) ax.set_xlabel(xaxis_label; fontsize = xaxis_label_fontsize); end;
    if !isnothing(yaxis_label) ax.set_ylabel(yaxis_label; fontsize = yaxis_label_fontsize); end;
    if !isnothing(xaxis_lims) ax.set_xlim(xaxis_lims); end;
    if !isnothing(yaxis_lims) ax.set_ylim(yaxis_lims); end;
    if !isnothing(xaxis_ticks) ax.set_xticks(xaxis_ticks); end;
    if !isnothing(yaxis_ticks) ax.set_yticks(yaxis_ticks); end;
    if !isnothing(xaxis_ticklabels) ax.set_xticklabels(xaxis_ticklabels); end;
    if !isnothing(yaxis_ticklabels) ax.set_yticklabels(yaxis_ticklabels); end;

    return nothing
);
