to_int(x::String) = parse(Int, x);

to_bool(x) = (
    if uppercase(x) in ["Y", "YES"]
        return true
    elseif uppercase(x) in ["N", "NO"]
        return false
    else
        return uppercase(x)
    end;
);
