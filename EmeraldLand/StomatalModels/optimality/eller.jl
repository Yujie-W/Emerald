# This file contains function to compute the partial derivatives of Θ versus E for EllerSM optimality model

∂Θ∂E!(cache::SPACCache{FT}, sm::EllerSM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}) where {FT} = (
    e = flow_out(leaf);
    δe = e / 100;
    dedp1 = ∂E∂P(leaf, e; δe = δe);
    dedp2 = ∂E∂P(leaf, e; δe = -δe);
    dkde  = (dedp2 - dedp1) / δe;
    leaf.flux.auxil.∂Θ∂E .= dkde .* leaf.flux.auxil.a_n ./ dedp1;

    return nothing
);
