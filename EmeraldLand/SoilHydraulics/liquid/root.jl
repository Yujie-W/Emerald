# This file contains function to account for root water uptake from the soil

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-13: add a utility function to read root water sink
# Note
#     This function is the same as the flow_in as in PlantHydraulics.jl. To avoid circular dependency, we copy and rename it to root_sink the function here.
#
#######################################################################################################################################################################################################
"""

    root_sink(root::Root{FT}) where {FT}

Return root water update, given
- `root` `Root` type struct that may contain non- and steady state flow

"""
function root_sink end;

root_sink(root::Root{FT}) where {FT} = root_sink(root.xylem.auxil);

root_sink(x_aux::XylemHydraulicsAuxilNSS{FT}) where {FT} = x_aux.flow[1];

root_sink(x_aux::XylemHydraulicsAuxilSS{FT}) where {FT} = x_aux.flow;


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Jun-29: tease apart this function for better readability
#     2023-Jul-06: add info into DEBUG code block
#
#######################################################################################################################################################################################################
"""

    root_source_sink!(spac::BulkSPAC{FT}) where {FT}

Update the source/sink terms for the soil layers, given
- `spac` the SPAC model

"""
function root_source_sink!(spac::BulkSPAC{FT}) where {FT}
    rindx = spac.plant.roots_index;
    roots = spac.plant.roots;
    sbulk = spac.soil_bulk;
    soils = spac.soils;

    # loop through the roots and compute the source/sink terms
    for i in eachindex(roots)
        soils[rindx[i]].auxil.∂θ∂t -= root_sink(roots[i]) * M_H₂O(FT) / ρ_H₂O(FT) / sbulk.trait.area / soils[rindx[i]].t_aux.δz;
    end;

    return nothing
end;
