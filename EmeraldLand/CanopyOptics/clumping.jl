#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-02: migrate the function from CanopyLayers
#     2022-Jun-02: rename the function from clumping_factor! to clumping_index!
#
#######################################################################################################################################################################################################
"""

    clumping_index!(can::MultiLayerCanopy) where {FT}

Update the clumping index, given
- `can` `MultiLayerCanopy` type canopy

"""
function clumping_index!(can::MultiLayerCanopy{FT}) where {FT}
    (; Ω_A, Ω_B) = can;

    can.ci = Ω_A + Ω_B * (1 - cosd(can.sun_geometry.state.sza));

    return nothing
end;
