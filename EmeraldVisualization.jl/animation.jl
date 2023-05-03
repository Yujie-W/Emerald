#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jan-20: add new function based on Plots.jl
#
#######################################################################################################################################################################################################
"""

    animation(make_frame!::Function, inds::Union{UnitRange{Int}, Vector{Int}}; filename::Union{Nothing,String} = nothing, fps::Int = 15)
    animation(dir::String, frames::Vector{String}; filename::Union{Nothing,String} = nothing, fps::Int = 15)

Generate animation, given
- `make_frame!` Function to make a frame
- `inds` Indices to pass to `make_frame!` function
- `filename` File name or path to save the animation (to gif or mp4; default is nothing)
- `fps` Frame per second (default is 15)
- `dir` Directory that saves the frames of an animation
- `frames` Frames of the animation (best to use PNG)

"""
function animation end

animation(make_frame!::Function, inds::Union{UnitRange{Int}, Vector{Int}}; filename::Union{Nothing,String} = nothing, fps::Int = 15) = (
    _anim = @animate for _i in inds
        make_frame!(_i);
    end;

    return isnothing(filename) ? gif(_anim; fps = fps) : (gif(_anim, filename; fps = fps); nothing)
);


animation(dir::String, frames::Vector{String}; filename::Union{Nothing,String} = nothing, fps::Int = 15) = (
    _anim = Animation(dir, frames);

    return isnothing(filename) ? gif(_anim; fps = fps) : (gif(_anim, filename; fps = fps); nothing)
);
