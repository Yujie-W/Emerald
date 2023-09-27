#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-25: migrate inititialize_legacy! function to clear_legacy!
#     2022-May-26: add method for each hydraulic system like leaf, root, and stem
#     2022-May-26: add method for leaf (hydraulic system nested within)
#     2022-May-26: add method for SPAC system (hydraulic system nested within)
#     2022-Jun-29: rename SPAC to ML*SPAC to be more accurate
#     2022-Jun-30: add compatibility to Leaves2D
#     2022-Jul-08: deflate documentations
#
#######################################################################################################################################################################################################
"""

    clear_legacy!(spac::MultiLayerSPAC{FT}) where {FT}
    clear_legacy!(organ::Union{Leaf2{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}) where {FT}

Clear the legacy for hydraulic organ or system, given
- `spac` `MultiLayerSPAC` type structure
- `organ` `Leaf2`, `Leaves2D`, `Root`, or `Stem` type structure
"""
function clear_legacy! end

clear_legacy!(spac::MultiLayerSPAC{FT}) where {FT} = (clear_legacy!.(spac.ROOTS); clear_legacy!(spac.TRUNK); clear_legacy!.(spac.BRANCHES); clear_legacy!.(spac.LEAVES););

clear_legacy!(organ::Union{Leaf2{FT}, Leaves2D{FT}, Root2{FT}, Stem2{FT}}) where {FT} = clear_legacy!(organ.NS);

clear_legacy!(xylem::XylemHydraulics{FT}) where {FT} = (xylem.auxil.k_history .= 1; xylem.state.p_history .= 0; return nothing);
