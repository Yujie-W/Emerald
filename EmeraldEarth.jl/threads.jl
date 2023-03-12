#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to load threadings
#
#######################################################################################################################################################################################################
"""

    add_threads!(nTH::Int)

Add processors to run code in multiple threadings, given
- `threads` Number of threads

"""
function add_threads!(threads::Int = 12)
    dynamic_workers!(threads);
    @everywhere Base.MainInclude.eval(using Emerald.EmeraldEarth);

    return nothing
end
