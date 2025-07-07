# This file contains function to compute the partial derivatives of Θ versus E for EllerSM optimality model
#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2023-Oct-16: make sure maximum gsc does not exceed g_CO₂_b
#     2024-Oct-16: make sure gsm is positive
#     2024-Oct-30: add leaf connection check
#
#######################################################################################################################################################################################################
∂Θ∂E!(cache::SPACCache{FT}, sm::EllerSM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}) where {FT} = (
    # if leaf xylem is not connected, do nothing
    if !leaf.xylem.state.connected
        leaf.flux.auxil.∂Θ∂E .= 0;

        return nothing
    end;

    # compute the ∂Θ∂E when leaf xylem is connected
    e = flow_out(leaf);
    δe = e / 100;
    dedp1 = ∂E∂P(leaf, e; δe = δe);
    dedp2 = ∂E∂P(leaf, e; δe = -δe);
    dkde  = (dedp2 - dedp1) / δe;
    leaf.flux.auxil.∂Θ∂E .= dkde .* leaf.flux.auxil.a_n ./ dedp1;

    return nothing
);
