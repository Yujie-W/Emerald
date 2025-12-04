"""
Struct that contains leaf traits and associated auxiliary vairables
"""
Base.@kwdef mutable struct LeafBio{FT<:AbstractFloat}
    # state variables (prognostic or structural)
    "Trait variables"
    trait::LeafBioTrait{FT} = LeafBioTrait{FT}()
    "State variables"
    state::LeafBioState{FT} = LeafBioState{FT}()
    "Auxiliary variables"
    auxil::LeafBioAuxil{FT}
end;

LeafBio(config::SPACConfig{FT}) where {FT} = return LeafBio{FT}(auxil = LeafBioAuxil(config));
