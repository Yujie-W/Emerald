#######################################################################################################################################################################################################
#
# Changes to this strcut
# General
#     2024-02-05: add struct for generalized flux tower data
#     2024-02-05: add constructor for default US_NR1 tower
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Generalized flux tower

$(TYPEDFIELDS)

"""
Base.@kwdef struct FluxTowerDataset
    "Flux tower label"
    LABEL::String
    "Temporal resolution"
    TRESO::String
    "Version tag of data"
    VER_TAG::String
    "Years of observation"
    YEARS::Vector{Int}
end;

FluxTowerDataset(tower::US_NR1) = FluxTowerDataset(
    LABEL = "AMF_US-NR1_BASE",
    TRESO = "HH",
    VER_TAG = "18-5",
    YEARS = collect(1998:2021),
);
