#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to load threadings
#
#######################################################################################################################################################################################################
"""

    add_threads!(threads::Int, FT::DataType)

Add processors to run code in multiple threadings, given
- `threads` Number of threads
- `FT` Floating type for the CACHE_SPAC

"""
function add_threads!(threads::Int, FT::DataType)
    @tinfo "Adding $(threads) threadings...";
    dynamic_workers!(threads);
    @everywhere Base.MainInclude.eval(:(using Emerald.EmeraldEarth));

    @tinfo "Initializing the SPAC cache in each thread...";
    @everywhere EmeraldEarth.initialize_cache!($FT);

    return nothing
end;
