#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-25: migrate inititialize_legacy! function to clear_legacy!
#     2022-May-26: add method for each hydraulic system like leaf, root, and stem
#     2022-May-26: add method for leaf (hydraulic system nested within)
#     2022-May-26: add method for SPAC system (hydraulic system nested within)
#     2022-Jun-29: rename SPAC to ML*SPAC to be more accurate
#     2022-Jun-30: add compatibility to Leaves1D and Leaves2D
#     2022-Jul-08: deflate documentations
#
#######################################################################################################################################################################################################
"""

    clear_legacy!(spac::MonoElementSPAC{FT}) where {FT}
    clear_legacy!(spac::MultiLayerSPAC{FT}) where {FT}
    clear_legacy!(organ::Union{Leaf{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}) where {FT}
    clear_legacy!(organ::Leaves1D{FT}) where {FT}

Clear the legacy for hydraulic organ or system, given
- `spac` `MonoElementSPAC` or `MultiLayerSPAC` type structure
- `organ` `Leaf`, `Leaves1D`, `Leaves2D`, `Root`, or `Stem` type structure
"""
function clear_legacy! end

clear_legacy!(spac::MonoElementSPAC{FT}) where {FT} = (clear_legacy!(spac.LEAF); clear_legacy!(spac.ROOT); clear_legacy!(spac.STEM););

clear_legacy!(spac::MultiLayerSPAC{FT}) where {FT} = (clear_legacy!.(spac.ROOTS); clear_legacy!(spac.TRUNK); clear_legacy!.(spac.BRANCHES); clear_legacy!.(spac.LEAVES););

clear_legacy!(organ::Union{Leaf{FT}, Leaves2D{FT}, Root{FT}, Stem{FT}}) where {FT} = clear_legacy!(organ.HS);

clear_legacy!(organ::Leaves1D{FT}) where {FT} = (clear_legacy!(organ.HS); clear_legacy!(organ.HS2););

clear_legacy!(hs::Union{LeafHydraulics{FT,DIM_XYLEM}, RootHydraulics{FT,DIM_XYLEM}, StemHydraulics{FT,DIM_XYLEM}}) where {FT,DIM_XYLEM} = (hs._k_history .= 1; hs.p_history .= 0; return nothing);
