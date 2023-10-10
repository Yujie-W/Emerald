# Operation functions for terminal IO
to_bool_or_uppercase(x) = (
    if uppercase(x) in ["Y", "YES"]
        return true
    elseif uppercase(x) in ["N", "NO"]
        return false
    else
        return uppercase(x)
    end;
);
to_nothing_or_int_vector(x::String) = (
    if x == ""
        return nothing
    elseif occursin(":", x)
        _x_split = split(x, ":");
        if length(_x_split) == 2
            _min = parse(Int, _x_split[1]);
            _max = parse(Int, _x_split[2]);
            return collect(_min:_max)
        else length(_x_split) > 2
            _min = parse(Int, _x_split[1]);
            _stp = parse(Int, _x_split[2]);
            _max = parse(Int, _x_split[3]);
            return collect(_min:_stp:_max)
        end;
    elseif occursin(",", x)
        _x_split = split(x, ",");
        return [parse(Int, _str) for _str in _x_split];
    else
        return [parse(Int, x)];
    end;
);
to_nothing_or_uppercase(x::String) = (x == "" ? nothing : uppercase(x));


# Judge functions for terminal IO
is_gm_mt(x::String) = has_no_space(x) && x[end;:end;] in ["H", "D", "M", "Y"];

is_nothing_or_all_letters(x::Union{Nothing,String}) = isnothing(x) || has_no_space(x);

is_nothing_or_vector(x) = isnothing(x) || x isa Vector;
