# This file contains function to test what happens if the xylem hydraulic system is in a new growth manner (not realistic, but has to be done this before I use segmented xylem VC)

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2024-Aug-30: add function to grow new xylem (area added to both total area and sap area)
#
#######################################################################################################################################################################################################
"""

    xylem_growth!(organ::Union{Root{FT}, Stem{FT}}, c_mol::FT)

Grow new xylem (area added to both total area and sap area), given
- `organ` `Root` or `Stem` type structure
- `c_mol` Carbon investment `[mol]`

"""
function xylem_growth! end;

xylem_growth!(organ::Union{Root{FT}, Stem{FT}}, c_mol::FT) where {FT} = xylem_growth!(organ.xylem, c_mol);

xylem_growth!(xylem::XylemHydraulics{FT}, c_mol::FT) where {FT} = (
    # compute the potential area from the given carbon investment
    delta_a = c_mol / (xylem.trait.l * xylem.trait.ρ * 1000 / 30);

    # add the new xylem to the xylem area and the sap area
    xylem.trait.area += delta_a;
    xylem.state.asap += delta_a;

    return nothing
);
