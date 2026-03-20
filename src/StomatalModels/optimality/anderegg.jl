"""

    ∂Θ∂E!(config::SPACConfig{FT}, cache::SPACCache{FT}, sm::AndereggSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT}

Update the marginal risk for stomatal opening, given
- `sm` `AndereggSM`, `EllerSM`, `SperrySM`, `Sperry2SM`, `WangSM`, or `Wang2SM` type optimality model
- `leaf` `Leaf` type struct

"""
∂Θ∂E!(config::SPACConfig{FT}, cache::SPACCache{FT}, sm::AndereggSM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
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
