module Math


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Oct-17: move function outside of the folder
#     2023-Jan-19: move function to EmeraldMath
#     2023-Mar-10: move function to Emerald
#
#######################################################################################################################################################################################################
"""

    numerical∫(f::Vector{FT}, Δx::Vector{FT}) where {FT}
    numerical∫(f::Vector{FT}, Δx::FT) where {FT}

Return the intergal of given
- `f` f(x) for each x
- `Δx` Δx for x

    numerical∫(f::Function, x_min::FT, x_max::FT, n::Int) where {FT}
    numerical∫(f::Function, x_min::FT, x_max::FT, x_tol::FT = sqrt(eps(FT)), y_tol::FT = sqrt(eps(FT))) where {FT}

Return the integral of given
- `f` A function
- `x_min` Minimum limit of x
- `x_max` Maximum limit of x
- `n` Number of points in the x range (evenly stepped)
- `x_tol` Tolerance of Δx (x/N)
- `y_tol` Tolerance of the integral solution

"""
function numerical∫ end

numerical∫(f::Vector{FT}, Δx::Vector{FT}) where {FT} = (
    if length(Δx) == length(f)
        return f' * Δx
    end;

    @warn "Dimensions not matching, use the matching parts only...";
    N = min(length(f), length(Δx));
    return view(f,1:N)' * view(Δx,1:N)
);

numerical∫(f::Vector{FT}, Δx::FT) where {FT} = sum(f) * abs(Δx);

numerical∫(f::Function, x_min::FT, x_max::FT, n::Int) where {FT} = (
    _sum = 0;
    _dx  = (x_max - x_min) / n;
    for _i in 1:n
        _x = x_min + FT(_i-0.5) * _dx;
        _sum += f(_x);
    end;

    return _sum * _dx
);

numerical∫(f::Function, x_min::FT, x_max::FT, x_tol::FT = sqrt(eps(FT)), y_tol::FT = sqrt(eps(FT))) where {FT} = (
    @assert y_tol > 0;

    # _sum_0: sum before halfing steps (_N), _sum_N: sum after halfing steps
    _sum_0::FT = (f(x_min) + f(x_max)) / 2;
    _sum_N::FT = 0;
    _N = 1;

    # continue the steps till the tolerances are reached
    while true
        _sum_N = (numerical∫(f, x_min, x_max, _N) + _sum_0) / 2;
        if abs(_sum_N - _sum_0) < y_tol || 1/_N < x_tol
            break;
        end;
        _sum_0 = _sum_N;
        _N *= 2;
    end;

    return _sum_N
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Oct-17: move functions outside of the folder
#     2023-Jan-19: move functions to EmeraldMath
#     2023-Mar-10: move function to Emerald
#
#######################################################################################################################################################################################################
"""

    lower_quadratic(a::FT, b::FT, c::FT) where {FT}

Return the lower quadratic solution or NaN, given
- `a` Parameter in `a*x^2 + b*x + c = 0`
- `b` Parameter in `a*x^2 + b*x + c = 0`
- `c` Parameter in `a*x^2 + b*x + c = 0`

"""
function lower_quadratic(a::FT, b::FT, c::FT) where {FT}
    discr = b^2 - 4*a*c;

    if (discr >= 0) && (a > 0)
        return (-b - sqrt(discr)) / (2 * a)
    elseif (discr >= 0) && (a < 0)
        return (-b + sqrt(discr)) / (2 * a)
    else
        return FT(NaN)
    end
end


"""

    upper_quadratic(a::FT, b::FT, c::FT) where {FT}

Return the upper quadratic solution or NaN, given
- `a` Parameter in `a*x^2 + b*x + c = 0`
- `b` Parameter in `a*x^2 + b*x + c = 0`
- `c` Parameter in `a*x^2 + b*x + c = 0`

"""
function upper_quadratic(a::FT, b::FT, c::FT) where {FT}
    discr = b^2 - 4*a*c;

    if (discr >= 0) && (a > 0)
        return (-b + sqrt(discr)) / (2 * a)
    elseif (discr >= 0) && (a < 0)
        return (-b - sqrt(discr)) / (2 * a)
    else
        return FT(NaN)
    end
end


end # Math
