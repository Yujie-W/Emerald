module StructEqual


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-27: add function compare_struct! to compare two structs
#     2024-Feb-28: add option show_diff_msg for testing purpose
#
#######################################################################################################################################################################################################
"""

    compare_struct!(struct1::ST, struct2::ST; first_element_array::Bool = true, show_diff_msg::Bool = true) where ST
    compare_struct!(struct1::ST, struct2::ST, target::Symbol; first_element_array::Bool = true, show_diff_msg::Bool = true) where ST

Comparing two structs `struct1` and `struct2` recursively, given
- `struct1` The first struct to compare
- `struct2` The second struct to compare
- `first_element_array` Whether to compare the first element of an array only
- `show_diff_msg` Whether to show the difference message
- `target` The field to compare

"""
function compare_struct! end;

compare_struct!(struct1::ST, struct2::ST; first_element_array::Bool = false, show_diff_msg::Bool = true) where ST = (
    n_error::Int = 0;

    for fn in fieldnames(ST)
        fntype = fieldtype(ST, fn);
        if fntype <: Union{Number, Bool}
            if !(getfield(struct1, fn) ≈ getfield(struct2, fn)) && !(isnan(getfield(struct1, fn)) && isnan(getfield(struct2, fn)))
                if show_diff_msg println("Field ", fn, " of ", ST, " is different!"); end;
                n_error += 1;
            end;
        elseif fntype <: String
            if !(getfield(struct1, fn) == getfield(struct2, fn))
                if show_diff_msg println("Field ", fn, " of ", ST, " is different!"); end;
                n_error += 1;
            end;
        elseif fntype <:AbstractArray && eltype(fntype) <: Union{Number, Bool}
            if !(getfield(struct1, fn) ≈ getfield(struct2, fn))
                if show_diff_msg println("Field ", fn, " of ", ST, " is different!"); end;
                n_error += 1;
            end;
        elseif fntype <:AbstractArray
            if first_element_array
                n_error += compare_struct!(getfield(struct1, fn)[1], getfield(struct2, fn)[1]; first_element_array = first_element_array, show_diff_msg = show_diff_msg);
            else
                for i in 1:length(getfield(struct1, fn))
                    n_error += compare_struct!(getfield(struct1, fn)[i], getfield(struct2, fn)[i]; first_element_array = first_element_array, show_diff_msg = show_diff_msg);
                end;
            end;
        elseif fntype <: Function
            nothing;
        else
            n_error += compare_struct!(getfield(struct1, fn), getfield(struct2, fn); first_element_array = first_element_array, show_diff_msg = show_diff_msg);
        end;
    end;

    return n_error
);

compare_struct!(struct1::ST, struct2::ST, target::Symbol; first_element_array::Bool = false, show_diff_msg::Bool = true) where ST = (
    n_error::Int = 0;

    if ST <: Union{Number, String, Bool, Function}
        return nothing
    end;

    if ST <: AbstractArray
        if eltype(ST) <: Union{Number, String, Bool, Function}
            return nothing
        else
            if first_element_array
                n_error += compare_struct!(struct1[1], struct2[1], target; first_element_array = first_element_array, show_diff_msg = show_diff_msg);
            else
                for i in 1:length(getfield(struct1, target))
                    n_error += compare_struct!(struct1[i], struct2[i], target; first_element_array = first_element_array, show_diff_msg = show_diff_msg);
                end;
            end;
        end;
    end;

    for fn in fieldnames(ST)
        if fn == target
            if show_diff_msg println("Comparing ", fn, " of ", ST); end;
            n_error += compare_struct!(getfield(struct1, fn), getfield(struct2, fn); first_element_array = first_element_array, show_diff_msg = show_diff_msg);
        else
            n_error += compare_struct!(getfield(struct1, fn), getfield(struct2, fn), target; first_element_array = first_element_array, show_diff_msg = show_diff_msg);
        end;
    end;

    return n_error
);


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
    for fn in fieldnames(ST)
        fntype = fieldtype(ST, fn);
        # TODO: memory allocation when sync numbers and bools
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
