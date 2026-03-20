∂Θ∂E!(config::SPACConfig{FT}, cache::SPACCache{FT}, sm::Sperry2SM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
    # if leaf xylem is not connected, do nothing
    if !leaf.xylem.state.connected
        leaf.flux.auxil.∂Θ∂E .= 0;

        return nothing
    end;

    (; NEW_C4_STOMATAL_METHODS) = config.FEATURES;

    # compute the ∂Θ∂E when leaf xylem is connected
    e = flow_out(leaf);
    δe = e / 100;
    dedp1 = ∂E∂P(leaf, e; δe = δe);
    dedp2 = ∂E∂P(leaf, e; δe = -δe);
    dedpm = ∂E∂P(leaf, FT(0); δe = δe);
    dkde  = (dedp2 - dedp1) / δe;
    am = NEW_C4_STOMATAL_METHODS ? min.(leaf.photosystem.auxil.a_p, leaf.photosystem.auxil.a_j) .- leaf.photosystem.auxil.r_d : leaf.flux.auxil.a_n;

    leaf.flux.auxil.∂Θ∂E .= dkde .* max.(FT(0.01), am) ./ dedpm;

    return nothing
);
