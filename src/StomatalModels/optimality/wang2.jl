# This file contains functions to compute the ∂Θ∂E for the Wang2SM optimization model (a modification of the AndereggSM optimization model)
#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-11: add method for Wang2SM model
#     2023-Mar-02: add a eps(FT) controller to (e_crit - e)
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2024-Oct-30: add leaf connection check
#
#######################################################################################################################################################################################################
∂Θ∂E!(cache::SPACCache{FT}, sm::Wang2SM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}) where {FT} = (
    # if leaf xylem is not connected, do nothing
    if !leaf.xylem.state.connected
        leaf.flux.auxil.∂Θ∂E .= 0;

        return nothing
    end;

    # compute the ∂Θ∂E when leaf xylem is connected
    dedp = ∂E∂P(leaf, flow_out(leaf)) / leaf.xylem.trait.area;
    leaf.flux.auxil.∂Θ∂E .= -sm.A .* leaf.capacitor.state.p_leaf .* leaf.flux.auxil.a_n ./ dedp;

    return nothing
);
