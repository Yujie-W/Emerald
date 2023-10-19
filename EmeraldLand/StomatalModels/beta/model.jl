# This file contains the functions to compute beta factor for empirical stomatal conductance models

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-30: migrate function from older version StomataModels.jl
#     2022-Jul-01: add method to tune stomatal opening based on relative hydraulic conductance at leaf xylem end;
#     2022-Jul-01: add method to tune stomatal opening based on relative hydraulic conductance of the soil
#     2022-Jul-01: add method to tune stomatal opening based on soil potential or leaf pressure
#     2022-Nov-18: force the beta to be within (0,1]
#
#######################################################################################################################################################################################################
"""

    β_factor(f::Function, x::FT) where {FT}

Return the β factor based on relative conductance or soil potential/pressure, given
- `f` Function to translate relative k to β, for example f(x) = x, f(x) = x², and f(x) = sqrt(x) for x in [0,1]
- `x` X variable to compute the β factor

Here, we pass a `f` (Function) into the function call so that we can easily customized the beta function without hacking the source code.

"""
function β_factor(f::Function, x::FT) where {FT}
    return max(eps(FT), min(1, f(x)))
end;
