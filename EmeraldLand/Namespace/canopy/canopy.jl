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


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2024-Feb-22: add struct MultiLayerCanopyStates
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save multiple layer hyperspectral canopy states (collections of states)

# Fields

$(TYPEDFIELDS)

"""
mutable struct MultiLayerCanopyStates{FT<:AbstractFloat}
    "Sensor geometry state"
    sensor_geometry::SensorGeometryState{FT}
    "Sun geometry state"
    sun_geometry::SunGeometryState{FT}
end;

MultiLayerCanopyStates(canopy::MultiLayerCanopy{FT}) where {FT} = MultiLayerCanopyStates{FT}(canopy.sensor_geometry.state, canopy.sun_geometry.state);

sync_state!(canopy::MultiLayerCanopy{FT}, states::MultiLayerCanopyStates{FT}) where {FT} = (
    sync_state!(canopy.sensor_geometry.state, states.sensor_geometry);
    sync_state!(canopy.sun_geometry.state, states.sun_geometry);

    return nothing
);

sync_state!(states::MultiLayerCanopyStates{FT}, canopy::MultiLayerCanopy{FT}) where {FT} = (
    sync_state!(states.sensor_geometry, canopy.sensor_geometry.state);
    sync_state!(states.sun_geometry, canopy.sun_geometry.state);

    return nothing
);
