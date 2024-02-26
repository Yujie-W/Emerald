# general function to update auxilary variables based on trait and state variables

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-25: add function t_aux!
#
#######################################################################################################################################################################################################
"""

    t_aux!(trait, t_aux)

Update the trait-dependent auxilary variables, given
- `trait` the trait variables
- `t_aux` the t_aux variables to be updated

"""
function t_aux! end;

t_aux!(config::SPACConfiguration, trait::Union{DataType,Nothing}, t_aux::Nothing) = nothing;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-25: add function s_aux!
#
#######################################################################################################################################################################################################
"""

    s_aux!(state, s_aux)

Update the state-dependent auxilary variables, given
- `state` the state variables
- `s_aux` the s_aux variables to be updated

"""
function s_aux! end;

s_aux!(config::SPACConfiguration, state::Union{DataType,Nothing}, s_aux::Nothing) = nothing;
