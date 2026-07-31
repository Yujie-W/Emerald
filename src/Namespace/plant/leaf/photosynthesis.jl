LeafPhotosystem(config::SPACConfig{FT}, c3c4::String = "C3") where {FT} = (
    @assert c3c4 in ["C3", "C4"] "The model string should be either C3 or C4!";

    cache_dim_ppar = isnothing(config.DIMENSIONS.DIM_PPAR_BINS) ? config.DIMENSIONS.DIM_INCL * config.DIMENSIONS.DIM_AZI : config.DIMENSIONS.DIM_PPAR_BINS;

    return if c3c4 == "C3"
        LeafPhotosystem{FT}(trait = C3Trait{FT}(), state = C3State{FT}(), auxil = LeafPhotosystemAuxil{FT}(cache_dim_ppar+1));
    else
        ps = LeafPhotosystem{FT}(trait = C4Trait{FT}(), state = C4State{FT}(), auxil = LeafPhotosystemAuxil{FT}(cache_dim_ppar+1));
        ps.auxil.f_psii = 0.41;
        ps
    end;
);
