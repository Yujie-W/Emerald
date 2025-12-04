#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the CACHE_SPAC
#     2023-Mar-13: initialize CACHE_CONFIG at the same time
#     2023-Jun-15: make sure prescribed swc does not exceed the limits
#     2024-Feb-22: remove state and auxil from spac struct
#     2024-Feb-23: rename to setup_cache!
#     2024-Feb-28: set up global config only
# Bug fixes
#     2023-Aug-26: make sure sza < 89 when total radiation is higher than 10 W m⁻²
#
#######################################################################################################################################################################################################
"""

    setup_cache!(FT::DataType = Float64)

Initialize the global parameter `CACHE_CONFIG` (in all threads after loading workers), given
- `FT` Floating type (default is Float64)

"""
function setup_cache!(FT::DataType = Float64)
    global CACHE_CONFIG = SPACConfiguration(FT);

    return nothing
end;
