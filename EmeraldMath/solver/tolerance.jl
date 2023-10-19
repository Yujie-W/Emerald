#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2023-Mar-10: move from ConstrainedRootSolvers to Emerald
#
#######################################################################################################################################################################################################
"""

Abstract tolerance type

$(TYPEDEF)

"""
abstract type AbstractTolerance{FT} end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-10: move from ConstrainedRootSolvers to Emerald
#
#######################################################################################################################################################################################################
"""

Tolerance for target function residual

$(TYPEDEF)

# Fields

$(TYPEDFIELDS)

"""
struct ResidualTolerance{FT} <: AbstractTolerance{FT}
    "Tolerance for residual"
    tol::FT
    "limit of iterations"
    n_limit::Int
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-10: move from ConstrainedRootSolvers to Emerald
#
#######################################################################################################################################################################################################
"""

Tolerance for solution

$(TYPEDEF)

# Fields

$(TYPEDFIELDS)

"""
struct SolutionTolerance{FT} <: AbstractTolerance{FT}
    "Tolerance for solution"
    tol::FT
    "limit of iterations"
    n_limit::Int
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-10: move from ConstrainedRootSolvers to Emerald
#
#######################################################################################################################################################################################################
"""

Tolerance for 2D and above solution

$(TYPEDEF)

# Fields

$(TYPEDFIELDS)

"""
struct SolutionToleranceND{FT} <: AbstractTolerance{FT}
    "Tolerance for solution"
    tol::Vector{FT}
    "limit of iterations"
    n_limit::Int
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: move function to Emerald
#
#######################################################################################################################################################################################################
"""

    if_break(tol::ResidualTolerance{FT}, x1::FT, x2::FT, y::FT, n::Int) where {FT}
    if_break(tol::SolutionTolerance{FT}, x1::FT, x2::FT, y::FT, n::Int) where {FT}

Determine whether to break, given
- `tol` tolerance struct
- `x1` Lower bound of x
- `x2` Upper bound of x
- `y` Residual of y
- `n` Current iteration

When the tolerance is for target function residual, if the modeled residual is lower than the given tolerance, or if the iteration exceeds the maximum limit, a `true` will be returned. When the
    tolerance is for solution, if the solution range is lower than the given tolerance, or if the iteration exceeds the maximum  limit, a `true` will be returned.

"""
function if_break end;

if_break(tol::ResidualTolerance{FT}, x1::FT, x2::FT, y::FT, n::Int) where {FT} = (abs(y) < tol.tol) || (n > tol.n_limit);

if_break(tol::SolutionTolerance{FT}, x1::FT, x2::FT, y::FT, n::Int) where {FT} = (abs(x1-x2) < tol.tol) || (n > tol.n_limit);
