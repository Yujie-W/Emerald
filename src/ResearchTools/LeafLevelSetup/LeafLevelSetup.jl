module LeafLevelSetup

using ..Namespace: SPACCache, SPACConfig
using ..Namespace: C3State, C4State, GeneralC3Trait, GeneralC4Trait, Leaf, LeafPhotosystem, LeafPhotosystemAuxil


include("cache.jl");
include("config.jl");
include("leaf.jl");
include("photosystem.jl");


end # module
