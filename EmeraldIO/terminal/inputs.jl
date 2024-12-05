#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-May-12: add function to read inputs from terminal
#     2023-Jul-14: add alias functions using the embedded operation and judge functions
#
#######################################################################################################################################################################################################
"""

    verified_input(message::String, operation_function::Function, judge_function::Function)
    verified_input(message::String, judge_function::Function)

Return a verified input, given
- `msg` Message about the input variable
- `operation_function` Function applied to the input before judge
- `judge_function` Function to judge if the input meets the requirements

"""
function verified_input end;

verified_input(message::String, operation_function::Function, judge_function::Function) = (
    _input = nothing;

    # loop until a correct input is given
    while true
        print(message);
        try
            _input = operation_function(readline());
            if judge_function(_input)
                break;
            else
                @warn "You input does not meet the requirements, please redo it!"
            end;
        catch e
            @warn e;
        end;
    end;

    return _input
);

verified_input(message::String, judge_function::Function) = (
    _input = nothing;

    # loop until a correct input is given
    while true
        print(message);
        _input = readline();
        if judge_function(_input)
            break;
        else
            @warn "You input does not meet the requirements, please redo it!"
        end;
    end;

    return _input
);

# alias functions
input_integer(msg::String; int_conversion::Bool = false, non_negative::Bool = false, positive::Bool = true) = (
    if positive
        _input = verified_input(msg, to_int, is_positive);
    elseif non_negative
        _input = verified_input(msg, to_int, is_not_negative);
    else
        _input = verified_input(msg, to_int);
    end;

    if int_conversion
        return parse(Int64, _input)
    else
        return _input
    end;
);

input_string(msg::String, operation_function::Function = (x -> x); allow_blank::Bool = false, no_space::Bool = false) = (
    if !allow_blank && no_space
        return verified_input(msg, operation_function, has_no_space)
    elseif !allow_blank && !no_space
        return verified_input(msg, operation_function, is_not_blank)
    elseif allow_blank && no_space
        _f(x) = has_no_space(x) || length(x) == 0;

        return verified_input(msg, operation_function, _f)
    else
        _g(x) = true;
        return verified_input(msg, operation_function, _g)
    end;
);

input_yes_or_no(msg::String; bool_conversion::Bool = false) = (
    _input = verified_input(msg, uppercase, is_yes_or_no);

    if bool_conversion
        return _input in ["Y", "YES"]
    else
        return _input
    end;
);
