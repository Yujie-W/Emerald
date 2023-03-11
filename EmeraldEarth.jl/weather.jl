#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to prescribe the weather drivers
#
#######################################################################################################################################################################################################
"""
"""
function prescribe! end

prescribe!(mat_spac::Matrix{Union{Nothing,MonoMLTreeSPAC}}, wd_tag::String; swc::Bool = true, t_leaf::Bool = true, t_soil::Bool = true) = (
    @assert wd_tag in ["wd1"] "Weather driver tag not supported: $(wd_tag)";

    # prescribe air layer environments

    return nothing
);
