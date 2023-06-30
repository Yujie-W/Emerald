#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-13: add a utility function to read root water sink
#
#######################################################################################################################################################################################################
"""

    root_sink(root::Root{FT}) where {FT<:AbstractFloat}

Return root water update, given
- `root` `Root` type struct that may contain non- and steady state flow

"""
function root_sink end

root_sink(root::Root{FT}) where {FT<:AbstractFloat} = root_sink(root.HS.FLOW);

root_sink(mode::SteadyStateFlow{FT}) where {FT<:AbstractFloat} = mode.flow;

root_sink(mode::NonSteadyStateFlow{FT}) where {FT<:AbstractFloat} = mode.f_in;
