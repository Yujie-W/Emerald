module StructEqual


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-24: add function sync_struct! to sync fields from one struct to another recursively
#
#######################################################################################################################################################################################################
"""

    sync_struct!(struct_from, struct_to)

Sync the fields from `struct_from` to `struct_to`.

"""
function sync_struct!(struct_from::ST, struct_to::ST) where ST
    for fn in fieldnames( ST )
        fntype = fieldtype(ST, fn);
        if fntype <: Union{Number, String, Bool}
            setfield!(struct_to, fn, getfield(struct_from, fn));
        elseif fntype <:AbstractArray && eltype(fntype) <: Union{Number, String, Bool}
            setfield!(struct_to, fn, getfield(struct_from, fn));
        elseif fntype <: Function
            setfield!(struct_to, fn, deepcopy(getfield(struct_from, fn)));
        else
            sync_struct!(getfield(struct_from, fn), getfield(struct_to, fn));
        end;
    end;

    return nothing
end;


end # module
