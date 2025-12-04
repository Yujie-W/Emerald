# This file contains functions to prescribe soil conditions

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-19: add function to update or prescribe swcs
#     2022-Oct-19: add option to prescribe t_soil profile
#     2023-Mar-28: update total energy in soil and leaf when prescribing swc and temperature
#     2023-Jun-12: update soil trace gas as well
#     2023-Jun-13: update N₂ and O₂ based on soil water content
#     2023-Jun-13: add soil gas energy into soil e
#     2023-Jun-15: make sure prescribed swc does not exceed the limits
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Oct-07: add 0.01 to the water vapor volume per soil layer
#     2023-Oct-18: re-initialize soil layer when any of the soil water content or temperature changes
#
#######################################################################################################################################################################################################

"""

    prescribe_soil!(spac::BulkSPAC{FT}; swcs::Union{Tuple,Nothing} = nothing, t_soils::Union{Tuple,Nothing} = nothing) where {FT}

Update the physiological parameters of the SPAC, given
- `config` Configuration for `BulkSPAC`
- `swcs` Soil water content at different layers. Optional, default is nothing
- `t_soils` Soil temperature at different layers. Optional, default is nothing

"""
function prescribe_soil!(spac::BulkSPAC{FT}; swcs::Union{Tuple,Nothing} = nothing, t_soils::Union{Tuple,Nothing} = nothing) where {FT}
    airs = spac.airs;
    soils = spac.soils;

    # prescribe soil water content (within [Θ_RES,Θ_SAT])
    if !isnothing(swcs)
        for i in eachindex(swcs)
            soils[i].state.θ = max(soils[i].trait.vc.Θ_RES + eps(FT), min(soils[i].trait.vc.Θ_SAT - eps(FT), swcs[i]));
        end;
    end;

    # prescribe soil temperature
    if !isnothing(t_soils)
        for i in eachindex(t_soils)
            soils[i].s_aux.t = t_soils[i];
        end;
    end;

    # if any of soil water content or temperature changes, re-initialize the soil layer
    if !isnothing(swcs) || !isnothing(t_soils)
        for soil in soils
            initialize_energy_states!(soil, airs[1]);
        end;
    end;

    return nothing
end;
