# general function to use in the namespace to synchronize the state variables between two structs

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-24: add function sync_state!
#
#######################################################################################################################################################################################################
"""

    sync_state!(from, to)

Synchonize the state variables from one to another

"""
function sync_state! end;

sync_state!(from::ST, to::ST) where {ST} = sync_struct!(from, to);
