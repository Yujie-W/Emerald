# This file contains functions to compute the partial derivatives of Θ versus E for AndereggSM optimality model

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: migrate function from older version
#     2022-Jul-11: add method for AndereggSM model
#     2022-Jul-11: add method for EllerSM model
#     2022-Jul-11: add method for SperrySM model
#     2022-Jul-11: add method for WangSM model
#     2022-Jul-11: add method for Wang2SM model
#     2023-Mar-02: add a eps(FT) controller to (e_crit - e)
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2024-Oct-30: add leaf connection check
#
#######################################################################################################################################################################################################
"""

    ∂Θ∂E!(cache::SPACCache{FT}, sm::AndereggSM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}) where {FT}

Update the marginal risk for stomatal opening, given
- `sm` `AndereggSM`, `EllerSM`, `SperrySM`, `WangSM`, or `Wang2SM` type optimality model
- `leaf` `Leaf` type struct

"""
∂Θ∂E!(cache::SPACCache{FT}, sm::AndereggSM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}) where {FT} = (
    # if leaf xylem is not connected, do nothing
    if !leaf.xylem.state.connected
        leaf.flux.auxil.∂Θ∂E .= 0;

        return nothing
    end;

    # compute the ∂Θ∂E when leaf xylem is connected
    dedp = ∂E∂P(leaf, flow_out(leaf)) / leaf.xylem.trait.area;
    leaf.flux.auxil.∂Θ∂E .= (-2 * sm.A * leaf.capacitor.state.p_leaf + sm.B) / dedp;

    return nothing
);
