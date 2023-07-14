has_no_space(x::String) = length(x) > 0 && !occursin(" ", x);

is_bool(x) = x isa Bool;

is_int(x::String) = (
    try
        parse(Int64, x);
        true;
    catch e
        false;
    end;
);

is_not_blank(x::String) = length(x) > 0 && any(!isspace, x);

is_not_negative(x::Number) = x >= 0;

is_positive(x::Number) = x > 0;

is_yes_or_no(x::String) = uppercase(x) in ["N", "NO", "Y", "YES"];
