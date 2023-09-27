#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-13: add a utility function to read root water sink
#
#######################################################################################################################################################################################################
"""

    root_sink(root::Root{FT}) where {FT}

Return root water update, given
- `root` `Root` type struct that may contain non- and steady state flow

"""
function root_sink end

root_sink(root::Root{FT}) where {FT} = root_sink(root.HS.FLOW);

root_sink(mode::SteadyStateFlow{FT}) where {FT} = mode.flow;

root_sink(mode::NonSteadyStateFlow{FT}) where {FT} = mode.f_in;


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Jun-29: tease apart this function for better readability
#     2023-Jul-06: add info into DEBUG code block
#
#######################################################################################################################################################################################################
"""

    soil_source_sink!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Update the source/sink terms for the soil layers, given
- `config` the SPAC configuration
- `spac` the SPAC model

"""
function soil_source_sink!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    (; DEBUG) = config;
    (; ROOTS, ROOTS_INDEX, SOIL) = spac;
    LAYERS = SOIL.LAYERS;

    # loop through the roots and compute the source/sink terms
    for _i in eachindex(ROOTS)
        LAYERS[ROOTS_INDEX[_i]].∂θ∂t -= root_sink(ROOTS[_i]) * M_H₂O(FT) / ρ_H₂O(FT) / SOIL.AREA / LAYERS[ROOTS_INDEX[_i]].ΔZ;
        LAYERS[ROOTS_INDEX[_i]].∂e∂t -= root_sink(ROOTS[_i]) / SOIL.AREA * CP_L_MOL(FT) * LAYERS[_i].t;

        if DEBUG
            if any(isnan, (LAYERS[ROOTS_INDEX[_i]].∂θ∂t, LAYERS[ROOTS_INDEX[_i]].∂e∂t))
                @info "Debugging" LAYERS[ROOTS_INDEX[_i]].∂θ∂t LAYERS[ROOTS_INDEX[_i]].∂e∂t;
                error("NaN detected in soil_source_sink! at layer $(ROOTS_INDEX[_i])");
            end;
        end;
    end;

    return nothing
end
