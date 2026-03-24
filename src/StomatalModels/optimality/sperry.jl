∂Θ∂E!(config::SPACConfig{FT}, cache::SPACCache{FT}, sm::SperrySM{FT}, leaf::Leaf{FT}, air::AirLayer{FT}) where {FT} = (
    # if leaf xylem is not connected, do nothing
    if !leaf.xylem.state.connected
        leaf.flux.auxil.∂Θ∂E .= 0;

        return nothing
    end;

    # compute the ∂Θ∂E when leaf xylem is connected
    e = flow_out(leaf);
    δe = max(min(e / 100, FT(1e-6)), FT(1e-7));
    dedp1 = ∂E∂P(leaf, e; δe = δe);
    dedp2 = ∂E∂P(leaf, e; δe = -δe);
    dedpm = ∂E∂P(leaf, FT(0); δe = δe);
    dkde  = (dedp2 - dedp1) / δe;

    # compute maximum A
    f_dif = relative_diffusive_coefficient(leaf.energy.auxil.t);
    g_max = leaf.flux.trait.g_limits[2] * f_dif;

    # use the cache vars
    gh1 = cache.cache_incl_azi_1_1;
    ghm = cache.cache_incl_azi_1_2;
    gsm = cache.cache_incl_azi_1_3;
    gcm = cache.cache_incl_azi_1_4;
    gh1 .= 1 ./ (1 ./ leaf.flux.state.g_H₂O_s .+ 1 ./ (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    ghm .= gh1 .* leaf.xylem.auxil.e_crit / e;
    gsm .= 1 ./ max.(eps(FT), 1 ./ ghm .- 1 ./ (FT(1.35) * leaf.flux.auxil.g_CO₂_b));
    gsm .= min.(gsm, g_max);
    gcm .= 1 ./ (FT(1.6) ./ gsm .+ 1 ./ leaf.flux.auxil.g_CO₂_b .+ 1 ./ leaf.flux.auxil.g_m);
    am = photosynthesis_only!(config, cache, leaf.photosystem, air, gcm, leaf.flux.auxil.ppar);

    leaf.flux.auxil.∂Θ∂E .= dkde .* max.(FT(0.01), am) ./ dedpm;

    return nothing
);
