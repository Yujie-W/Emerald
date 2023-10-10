#=

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Mar-28: add function to disconnect root from soil
#     2023-Mar-28: add function to disconnect others (by setting flows to 0)
#     2023-Sep-14: set g of leaves to 0 when root system is disconnected
#
#######################################################################################################################################################################################################
"""

    disconnect!(organ::Union{Leaves2D{FT},Stem{FT}}) where {FT}
    disconnect!(organ::Root{FT}) where {FT}

Disconnect root from soil (and set othes' flow to 0), given
- `organ` Root, stem, or leaf

"""
function disconnect! end;

disconnect!(organ::Leaves2D{FT}) where {FT} = (
    disconnect!(organ.HS, organ.HS.FLOW);

    organ.g_H₂O_s_shaded = 0;
    organ.g_H₂O_s_sunlit .= 0;
    organ._g_CO₂_shaded = 0;
    organ._g_CO₂_sunlit .= 0;

    return nothing
);

disconnect!(organ::Root2{FT}) where {FT} = (
    organ._isconnected = false;
    disconnect!(organ.HS, organ.HS.FLOW);

    return nothing
);

disconnect!(organ::Stem2{FT}) where {FT} = (
    disconnect!(organ.HS, organ.HS.FLOW);

    return nothing
);

disconnect!(hs::LeafHydraulics{FT}, mode::NonSteadyStateFlow{FT}) where {FT} = (
    # update the pressure
    hs.p_leaf = hs._p_storage;

    # update the flow
    mode.f_in = 0;
    mode.f_out = 0;
    mode._f_buffer .= 0;
    mode._f_element .= 0;
    mode._f_sum .= 0;

    return nothing
);

disconnect!(hs::Union{RootHydraulics{FT}, StemHydraulics{FT}}, mode::NonSteadyStateFlow{FT}) where {FT} = (
    # update the pressure
    hs._p_element .= hs._p_storage;

    # update the flow
    mode.f_in = 0;
    mode.f_out = 0;
    mode._f_buffer .= 0;
    mode._f_element .= 0;
    mode._f_sum .= 0;

    return nothing
);

disconnect!(hs::Union{LeafHydraulics{FT}, RootHydraulics{FT}, StemHydraulics{FT}}, mode::SteadyStateFlow{FT}) where {FT} = (
    mode.flow = 0;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-12: add function to disconnect roots based on soil water content
#
#######################################################################################################################################################################################################
"""

    disconnect_roots!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Disconnect roots from soil if soil water content is very low, given
- `config` SPAC configuration
- `spac` SPAC model

"""
function disconnect_roots!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    (; KR_THRESHOLD) = config;
    (; ROOTS, ROOTS_INDEX, SOIL) = spac;

    # very first step here: if soil is too dry, disconnect root from soil
    for i in eachindex(ROOTS)
        _root = ROOTS[i];
        soil = SOIL.LAYERS[ROOTS_INDEX[i]];
        _ψ_soil = soil_ψ_25(soil.VC, soil.θ) * relative_surface_tension(soil.t);
        _p_crit = xylem_pressure(_root.HS.VC, KR_THRESHOLD) * relative_surface_tension(_root.t);
        if _ψ_soil <= _p_crit
            disconnect!(_root);
        else
            _root._isconnected = true;
        end;
    end;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-12: add function to disconnect spac if no root is connected to soil any more
#
#######################################################################################################################################################################################################
"""

    disconnect_spac!(spac::MultiLayerSPAC{FT}) where {FT}

Disconnect spac if no root is connected to soil any more, given
- `spac` SPAC model

"""
function disconnect_spac!(spac::MultiLayerSPAC{FT}) where {FT}
    (; BRANCHES, LEAVES, ROOTS, TRUNK) = spac;

    # very first step here: if soil is too dry, disconnect root from soil
    _connected = 0;
    for _root in ROOTS
        _connected += _root._isconnected;
    end;

    # if all roots are disconnected, set all flows to 0
    # TODO: if leaf shedding is allowed, remember to update leaf area when roots are reconnected to the soil
    if _connected > 0
        if !spac._root_connection
            @info " Root system is now reconnected to soil!";
        end;
        spac._root_connection = true;
    else
        if spac._root_connection
            @info "Root system is now all disconnected from soil!"
            spac._root_connection = false;
            disconnect!(TRUNK);
            disconnect!.(BRANCHES);
            disconnect!.(LEAVES);
        end;
    end;

    return nothing
end;

=#
