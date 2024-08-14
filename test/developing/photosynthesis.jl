# This script is meant to test the new photosynthesis model

using Emerald
using Revise
using Test


FT = Float64;

# C3 models
begin "Photosynthesis model with GeneralC3Trait"
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    ps = EmeraldLand.Namespace.LeafPhotosystem{FT}("C3VJP");
    air = EmeraldLand.Namespace.AirLayer{FT}();
    p_i = FT(20);
    ppar = FT(1000);
    t_leaf = FT(298.15);
    EmeraldLand.Photosynthesis.photosynthesis!(config, ps, air, p_i, ppar, t_leaf);
end;


begin "Photosynthesis model with GeneralC3Trait + Cytochrome"
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    ps = EmeraldLand.Namespace.LeafPhotosystem{FT}("C3Cyto");
    air = EmeraldLand.Namespace.AirLayer{FT}();
    p_i = FT(20);
    ppar = FT(1000);
    t_leaf = FT(298.15);
    EmeraldLand.Photosynthesis.photosynthesis!(config, ps, air, p_i, ppar, t_leaf);
end;
