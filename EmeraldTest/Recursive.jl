module Recursive


#######################################################################################################################################################################################################
#
# Changes to these functions
# General
#     2022-Aug-24: move function outside of the folder
#     2022-Aug-24: simply the algorithm
#
#######################################################################################################################################################################################################
"""

    FT_test(para::Array, FT)
    FT_test(para::Number, FT)
    FT_test(para::Union{Function, Module, String, Symbol}, FT)
    FT_test(para::Any, FT)

Return true or false to determine if the FT is consistent, given
- `para` Parameter to run FT control
- `FT` Float type

If the testing variable is an array, the function will test if element type is float number:
- If true, the function tests if the element type is the same as given `FT`
- If false, the function tests each element recursively

The variable to test maybe a struct, but `FT_test` does not know the struct type name a priori. Thus, we try to read out the fields of the variable:
- If succeeds, the function test the fields recursively
- If fails, then do nothing

---
Example
```julia
struct SA
    a
    b
end;
sa = SA(1, 2.0);

ft_1 = FT_test([1, 2, 3], Float64);
ft_2 = FT_test(Any[1, 1.0f0, 1.0e0], Float64);
ft_3 = FT_test([1, 2.0, "a"], Float64);
ft_4 = FT_test(sa, Float64);
```

"""
function FT_test end;

FT_test(para::Array, FT) = (
    # fail if para is float but not FT
    if eltype(para) <: AbstractFloat
        return eltype(para) == FT
    end;

    # test all the elements
    return all(FT_test.(para, FT))
);

FT_test(para::Number, FT) = (
    # fail if para is float but not FT
    if typeof(para) <: AbstractFloat
        return typeof(para) == FT
    end;

    return true
);

FT_test(para::Union{Function, Module, String, Symbol}, FT) = true;

FT_test(para::Any, FT) = (
    # try to detech struct
    if !(typeof(para) <: DataType)
        try
            arr = [];
            for fn in fieldnames( typeof(para) )
                push!(arr, FT_test(getfield(para, fn), FT));
            end;

            return all(arr)
        catch e
            nothing
        end;
    end;

    return true
);


#######################################################################################################################################################################################################
#
# Changes to these functions
# General
#     2022-Aug-24: move function outside of the folder
#     2022-Aug-24: simply the algorithm
#
#######################################################################################################################################################################################################
"""
Like [`FT_test`](@ref), same logic is used to test if all the elements within the tested variable are not NaN:

    NaN_test(para::Array)
    NaN_test(para::Number)
    NaN_test(para::Union{Function, Module, String, Symbol})
    NaN_test(para::Any)

Test if the variable is not NaN, given
- `para` Parameter to test

---
Example
```julia
struct SA
    a
    b
end;

nan_1 = NaN_test(SA(1,2));
nan_2 = NaN_test(SA(1,NaN));
nan_3 = NaN_test([1,2,NaN]);
nan_4 = NaN_test([1,3,4]);
nan_5 = NaN_test([1,2,"a"]);
```

"""
function NaN_test end;

NaN_test(para::Array) = all(NaN_test.(para));

NaN_test(para::Number) = !isnan(para);

NaN_test(para::Union{Function, Module, String, Symbol}) = true;

NaN_test(para::Any) = (
    # try to detech struct
    if !(typeof(para) <: DataType)
        try
            arr = [];
            for fn in fieldnames( typeof(para) )
                push!(arr, NaN_test( getfield(para, fn) ));
            end;

            return all(arr)
        catch e
            nothing
        end;
    end;

    return true
);


end;
