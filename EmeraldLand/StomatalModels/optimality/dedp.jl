# This file contains function to compute ∂E∂P to use with optimality models

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jul-08: migrate function from StomataModels.jl
#     2022-Jul-08: add method for Leaf
#     2022-Jul-08: add option δe to be more general
#
#######################################################################################################################################################################################################
"""

    ∂E∂P(leaf::Leaf{FT}, flow::FT; δe::FT = FT(1e-7)) where {FT}

Return the marginal hydraulic conductance, given
- `leaf` `Leaf` type struct
- `flow` Flow rate through the leaf xylem `[mol s⁻¹]`
- `δe` Incremental flow rate, default is 1e-7

"""
function ∂E∂P(leaf::Leaf{FT}, flow::FT; δe::FT = FT(1e-7)) where {FT}
    @assert δe != 0 "δe must not be 0";

    p1 = xylem_end_pressure(leaf.xylem, flow, leaf.energy.s_aux.t);
    p2 = xylem_end_pressure(leaf.xylem, flow + δe, leaf.energy.s_aux.t);
    dedp = -δe / (p2 - p1);

    return dedp
end;
