#=
#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to load threadings
#
#######################################################################################################################################################################################################
"""

    add_threads!(threads::Int, FT::DataType = Float64)

Add processors to run code in multiple threadings, given
- `threads` Number of threads
- `FT` Floating type for the SPAC (default is Float64)

"""
function add_threads!(threads::Int, FT::DataType = Float64)
    pretty_display!("Adding a total of $(threads) threadings...", "tinfo_pre");
    dynamic_workers!(threads);
    # @everywhere Base.MainInclude.eval(:(using Emerald.EmeraldEarth));
    @everywhere eval(:(using Emerald.EmeraldEarth));

    pretty_display!("Initializing the SPAC cache in each thread...", "tinfo_end");
    @everywhere EmeraldEarth.setup_cache!($FT);

    return nothing
end;
=#
