module StructEqual


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-27: add function compare_struct! to compare two structs
#     2024-Feb-28: add option show_diff_msg for testing purpose
#     2024-Jul-25: add approximation option for comparing numbers
#
#######################################################################################################################################################################################################
"""

    compare_struct!(struct1::ST, struct2::ST; approximation::Bool = true, first_element_array::Bool = true, show_diff_msg::Bool = true) where ST
    compare_struct!(struct1::ST, struct2::ST, target::Symbol; first_element_array::Bool = true, show_diff_msg::Bool = true) where ST

Comparing two structs `struct1` and `struct2` recursively, given
- `struct1` The first struct to compare
- `struct2` The second struct to compare
- `approximation` Whether to compare the numbers approximately
- `first_element_array` Whether to compare the first element of an array only
- `show_diff_msg` Whether to show the difference message
- `target` The field to compare

"""
function compare_struct! end;

compare_struct!(struct1::ST, struct2::ST; approximation::Bool = true, first_element_array::Bool = false, show_diff_msg::Bool = true) where ST = (
    n_error::Int = 0;

    for fn in fieldnames(ST)
        fntype = fieldtype(ST, fn);
        if fntype <: Union{Number, Bool}
            if approximation
                if !(getfield(struct1, fn) ≈ getfield(struct2, fn)) && !(isnan(getfield(struct1, fn)) && isnan(getfield(struct2, fn)))
                    if show_diff_msg
                        println("Field ", fn, " of ", ST, " is different!");
                        println("    struct1: ", getfield(struct1, fn));
                        println("    struct2: ", getfield(struct2, fn));
                    end;
                    n_error += 1;
                end;
            else
                if !(getfield(struct1, fn) == getfield(struct2, fn)) && !(isnan(getfield(struct1, fn)) && isnan(getfield(struct2, fn)))
                    if show_diff_msg
                        println("Field ", fn, " of ", ST, " is different!");
                        println("    struct1: ", getfield(struct1, fn));
                        println("    struct2: ", getfield(struct2, fn));
                    end;
                    n_error += 1;
                end;
            end;
        elseif fntype <: String
            if !(getfield(struct1, fn) == getfield(struct2, fn))
                if show_diff_msg
                    println("Field ", fn, " of ", ST, " is different!");
                    println("    struct1: ", getfield(struct1, fn));
                    println("    struct2: ", getfield(struct2, fn));
                end;
                n_error += 1;
            end;
        elseif fntype <:AbstractArray && eltype(getfield(struct1, fn)) <: Union{Number, Bool}
            if approximation
                if !(getfield(struct1, fn) ≈ getfield(struct2, fn))
                    if show_diff_msg
                        println("Field ", fn, " of ", ST, " is different!");
                    end;
                    n_error += 1;
                end;
            else
                if !(getfield(struct1, fn) == getfield(struct2, fn))
                    if show_diff_msg
                        println("Field ", fn, " of ", ST, " is different!");
                    end;
                    n_error += 1;
                end;
            end;
        elseif fntype <:AbstractArray
            if first_element_array
                n_error += compare_struct!(getfield(struct1, fn)[1], getfield(struct2, fn)[1]; approximation = approximation, first_element_array = first_element_array, show_diff_msg = show_diff_msg);
            else
                for i in 1:length(getfield(struct1, fn))
                    n_error += compare_struct!(getfield(struct1, fn)[i], getfield(struct2, fn)[i]; approximation = approximation, first_element_array = first_element_array, show_diff_msg = show_diff_msg);
                end;
            end;
        elseif fntype <: Function
            nothing;
        else
            n_error += compare_struct!(getfield(struct1, fn), getfield(struct2, fn); approximation = approximation, first_element_array = first_element_array, show_diff_msg = show_diff_msg);
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
        if eltype(struct1) <: Union{Number, String, Bool, Function}
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
#     2025-Jun-03: add method to sync the struct to a dictionary (to write to file)
#     2025-Jun-04: add method to sync the struct from a dictionary (to read from file)
#
#######################################################################################################################################################################################################
"""

    sync_struct!(struct_from, struct_to)

Sync the fields from `struct_from` to `struct_to`.

"""
function sync_struct! end;

sync_struct!(struct_from::ST, struct_to::ST) where ST = (
    for fn in fieldnames(ST)
        fntype = fieldtype(ST, fn);
        # TODO: memory allocation when sync numbers and bools
        if fntype <: Union{Number, String, Bool}
            setfield!(struct_to, fn, getfield(struct_from, fn));
        elseif fntype <:AbstractArray && eltype(getfield(struct_from, fn)) <: Union{Number, String, Bool}
            setfield!(struct_to, fn, getfield(struct_from, fn));
        elseif fntype <: Function
            setfield!(struct_to, fn, deepcopy(getfield(struct_from, fn)));
        else
            sync_struct!(getfield(struct_from, fn), getfield(struct_to, fn));
        end;
    end;

    return nothing
);

sync_struct!(struct_from::ST) where ST = (
    dict = Dict{String,Any}();

    for fn in fieldnames(ST)
        fntype = fieldtype(ST, fn);

        if fntype <: Union{Number, String, Bool}
            dict[String(fn)] = getfield(struct_from, fn);
        elseif fntype <:AbstractArray && eltype(getfield(struct_from, fn)) <: Union{AbstractArray, Number, String, Bool}
            dict[String(fn)] = getfield(struct_from, fn);
        elseif fntype <:AbstractArray
            dict[String(fn)] = [
                sync_struct!(getfield(struct_from, fn)[i]) for i in eachindex(getfield(struct_from, fn))
            ];
        elseif fntype <: Function
            @warn "Function field $fn is not copied...";
        else
            dict[String(fn)] = sync_struct!(getfield(struct_from, fn));
        end;
    end;

    return dict
);

sync_struct!(dict_from::Dict, struct_to::ST) where ST = (
    fns = fieldnames(ST);

    for (k,v) in dict_from
        if Symbol(k) in fns
            fntype = fieldtype(ST, Symbol(k));

            if fntype <: Union{Number, String, Bool}
                setfield!(struct_to, Symbol(k), v);
            elseif fntype <:AbstractArray && eltype(v) <: Union{Number, String, Bool}
                setfield!(struct_to, Symbol(k), v);
            elseif fntype <:AbstractArray
                for idict in eachindex(v)
                    sync_struct!(v[idict], getfield(struct_to, Symbol(k))[idict]);
                end;
            elseif fntype <: Function
                @warn "Function field $(k) is not copied...";
            else
                sync_struct!(v, getfield(struct_to, Symbol(k)));
            end;
        else
            @warn "Field $(k) not existing in the target struct!";
        end;
    end;

    return nothing
);


end # module
