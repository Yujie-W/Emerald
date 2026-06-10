"""
Hierarchy of `AbstractCanopyRT`
- `CanopyRTSCOPE`
- `CanopyRTEmerald`

"""
abstract type AbstractCanopyRT end;

""" Canopy RT method based on the SCOPE model (Soil Canopy Observation of Photosynthesis and Energy fluxes) """
struct CanopyRTSCOPE <: AbstractCanopyRT end;

""" Canopy RT method of the 4SAIL model adapted from SCOPE """
struct CanopyRTEmerald <: AbstractCanopyRT end;


# Union alias
UnionCanopyRTMethod = Union{CanopyRTSCOPE, CanopyRTEmerald};
