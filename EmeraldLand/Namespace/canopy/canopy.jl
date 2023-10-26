# This file conmtains the struct for the canopy structure

#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2023-Oct-14: add struct MultiLayerCanopy
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save multiple layer hyperspectral canopy parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MultiLayerCanopy{FT<:AbstractFloat}
    # sun-sensor geometry related structs
    "Sensor geometry information"
    sensor_geometry::SensorGeometry{FT}
    "Sun geometry information"
    sun_geometry::SunGeometry{FT}
    "Canopy structure"
    structure::CanopyStructure{FT}
end;

MultiLayerCanopy(config::SPACConfiguration{FT}, n_layer::Int) where {FT} = (
    return MultiLayerCanopy{FT}(
                sensor_geometry = SensorGeometry(config, n_layer),
                sun_geometry    = SunGeometry(config, n_layer),
                structure       = CanopyStructure(config, n_layer),
    )
);
