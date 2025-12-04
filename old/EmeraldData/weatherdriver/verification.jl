#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Feb-29: add function to verify if the weather driver contains NaN (non breaking)
#
#######################################################################################################################################################################################################
"""

    dict_contains_nan(dict::Dict{String,Any})

Verify if the dictionary contains NaN, given
- `dict` Dictionary to verify

"""
function dict_contains_nan(dict::Dict{String,Any})
    n_error::Int = 0;

    # iterate through the keys and values
    for (_, value) in dict
        if typeof(value) <: Number
            if isnan(value)
                n_error += 1;
            end;
        elseif typeof(value) <: Array
            if any(isnan.(value))
                n_error += 1;
            end;
        end;
    end;

    return n_error > 0
end;
