"""

    leaf_level_photosystem(FT, c3c4::String)

Create and return a `LeafPhotosystem` struct for leaf-level simulations, given
- `FT` Floating-point type (e.g., `Float32`, `Float64`)
- `c3c4` String indicating whether to create a C3 or C4 photosystem

"""
function leaf_level_photosystem(FT, c3c4::String)
    @assert c3c4 in ["C3", "C4"] "The model string should be either C3 or C4!";

    return if c3c4 == "C3"
        LeafPhotosystem{FT}(trait = GeneralC3Trait{FT}(), state = C3State{FT}(), auxil = LeafPhotosystemAuxil{FT}(1));
    else
        ps = LeafPhotosystem{FT}(trait = GeneralC4Trait{FT}(), state = C4State{FT}(), auxil = LeafPhotosystemAuxil{FT}(1));
        ps.auxil.f_psii = 0.41;
        ps
    end;
end;
