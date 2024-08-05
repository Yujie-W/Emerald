# This file contains functions to compute the ∂Θ∂E for the Wang2SM optimization model (a modification of the AndereggSM optimization model)

∂Θ∂E!(cache::SPACCache{FT}, sm::Wang2SM{FT}, leaf::CanopyLayer{FT}, air::AirLayer{FT}) where {FT} = (
    dedp = ∂E∂P(leaf, flow_out(leaf)) / leaf.xylem.trait.area;
    leaf.flux.auxil.∂Θ∂E .= -sm.A .* leaf.capacitor.state.p_leaf .* leaf.flux.auxil.a_n ./ dedp;

    return nothing
);
