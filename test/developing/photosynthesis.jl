# This script is meant to test the new photosynthesis model

using Emerald
using Revise
using Test


FT = Float64;

# C3 models
begin "Photosynthesis model with GeneralC3Trait"
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    ps = EmeraldLand.Namespace.LeafPhotosystem{FT}();
    ps.trait = EmeraldLand.Namespace.GeneralC3Trait{FT}();
    air = EmeraldLand.Namespace.AirLayer{FT}();
    p_i = FT(20);
    ppar = FT(1000);
    t_leaf = FT(298.15);
    EmeraldLand.Photosynthesis.photosynthesis!(config, ps, air, p_i, ppar, t_leaf);
end;


begin "Photosynthesis model with GeneralC3Trait + Cytochrome"
    config = EmeraldLand.Namespace.SPACConfiguration(FT);
    ps = EmeraldLand.Namespace.LeafPhotosystem{FT}();
    ps.trait = EmeraldLand.Namespace.GeneralC3Trait{FT}();
    ps.trait.AJM = EmeraldLand.Namespace.AjMethodC3VqmaxPi();
    ps.trait.FLM = EmeraldLand.Namespace.CytochromeFluoscenceModel{FT}();
    air = EmeraldLand.Namespace.AirLayer{FT}();
    p_i = FT(20);
    ppar = FT(1000);
    t_leaf = FT(298.15);
    EmeraldLand.Photosynthesis.photosynthesis!(config, ps, air, p_i, ppar, t_leaf);
end;
