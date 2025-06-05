# This file contains functions to balance the air volume in the soil

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Jun-30: move function out of soil_budget!
#     2023-Jul-06: sort the order of gas volume balance and water volume balance
#     2023-Jul-06: add PRESCRIBE_AIR mode to avoid the errors due to mass balance in air
#     2023-Oct-07: limit the volume change from the source to 1/2 of the total dry air (soil and air)
#     2025-Jun-05: make soil total energy relative to triple temperature for phase change purposes
#     2025-Jun-05: account for ice volume in the calculation of air volume balance
#
#######################################################################################################################################################################################################
"""

    volume_balance!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Balance the air volume in the soil so that pressure is in equilibrium, given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC

"""
function volume_balance!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    (; PRESCRIBE_AIR) = config;
    air = spac.airs[1];
    soils = spac.soils;
    top_soil = soils[1];

    # balance the air volume among soil layers from lower to upper layers
    for i in length(soils)-1:-1:1
        # upper layer is soil_i and lower is soil_j
        soil_i = soils[i];
        soil_j = soils[i+1];

        # compute the air moles in the lower layer
        ndry_i = soil_i.state.ns[1] + soil_i.state.ns[2] + soil_i.state.ns[4] + soil_i.state.ns[5];
        ndry_j = soil_j.state.ns[1] + soil_j.state.ns[2] + soil_j.state.ns[4] + soil_j.state.ns[5];
        nmax_j = (air.state.p_air - saturation_vapor_pressure(soil_j.s_aux.t, soil_j.s_aux.ψ * 1000000)) * soil_j.t_aux.δz *
                 (soil_j.trait.vc.Θ_SAT - soil_j.state.θ - soil_j.state.θ_ice) / (GAS_R(FT) * soil_j.s_aux.t);

        # if ndry_j == nmax_j, no air needs to be transferred from/to the upper layer

        # if nmax_j > ndry_j and ndry_i > 0, air needs to be transferred from the upper layer
        if (nmax_j > ndry_j) && (ndry_i > 0)
            n_mass = min(nmax_j - ndry_j, ndry_i / 2);
            soil_j.state.ns[1] += n_mass * soil_i.state.ns[1] / ndry_i;
            soil_j.state.ns[2] += n_mass * soil_i.state.ns[2] / ndry_i;
            soil_j.state.ns[4] += n_mass * soil_i.state.ns[4] / ndry_i;
            soil_j.state.ns[5] += n_mass * soil_i.state.ns[5] / ndry_i;
            soil_j.state.Σe += n_mass * CP_D_MOL(FT) * (soil_i.s_aux.t - T₀(FT)) / soil_j.t_aux.δz;
            soil_i.state.ns[1] -= n_mass * soil_i.state.ns[1] / ndry_i;
            soil_i.state.ns[2] -= n_mass * soil_i.state.ns[2] / ndry_i;
            soil_i.state.ns[4] -= n_mass * soil_i.state.ns[4] / ndry_i;
            soil_i.state.ns[5] -= n_mass * soil_i.state.ns[5] / ndry_i;
            soil_i.state.Σe -= n_mass * CP_D_MOL(FT) * (soil_i.s_aux.t - T₀(FT)) / soil_i.t_aux.δz;
        end;

        # if nmax_j > ndry_j but upper layer is saturated, the lower layer will need to suck some water from upper layer to balance the air volume
        # compute the equilibrate mole of air from upper layer when it reaches its field capacity
        # TODO: not sure if ice cause any trouble here, as θ will be lower than Θ_SAT when there is ice, so it should be fine
        if (nmax_j > ndry_j) && (soil_i.state.θ >= soil_i.trait.vc.Θ_SAT)
            θ_fc_up = soil_θ(soil_i.trait.vc, -1 * ρg_MPa(FT) * soil_i.t_aux.δz / relative_surface_tension(soil_i.s_aux.t));
            v_mass = min((nmax_j - ndry_j) * GAS_R(FT) * soil_j.s_aux.t / air.state.p_air, (soil_i.trait.vc.Θ_SAT - θ_fc_up) * soil_i.t_aux.δz);
            soil_j.state.θ += v_mass / soil_j.t_aux.δz;
            soil_i.state.θ -= v_mass / soil_i.t_aux.δz;
            soil_j.state.Σe += v_mass * ρ_H₂O(FT) * CP_L(FT) * (soil_i.s_aux.t - T₀(FT)) / soil_j.t_aux.δz;
            soil_i.state.Σe -= v_mass * ρ_H₂O(FT) * CP_L(FT) * (soil_i.s_aux.t - T₀(FT)) / soil_i.t_aux.δz;
        end;

        # if nmax_j < ndry_j and ndry_j > 0, air needs to be transferred to the upper layer (does not matter whether upper layer is saturated or not)
        if (nmax_j < ndry_j) && (ndry_j > 0)
            n_mass = ndry_j - nmax_j;
            soil_j.state.ns[1] -= n_mass * soil_j.state.ns[1] / ndry_j;
            soil_j.state.ns[2] -= n_mass * soil_j.state.ns[2] / ndry_j;
            soil_j.state.ns[4] -= n_mass * soil_j.state.ns[4] / ndry_j;
            soil_j.state.ns[5] -= n_mass * soil_j.state.ns[5] / ndry_j;
            soil_j.state.Σe -= n_mass * CP_D_MOL(FT) * (soil_j.s_aux.t - T₀(FT)) / soil_j.t_aux.δz;
            soil_i.state.ns[1] += n_mass * soil_j.state.ns[1] / ndry_j;
            soil_i.state.ns[2] += n_mass * soil_j.state.ns[2] / ndry_j;
            soil_i.state.ns[4] += n_mass * soil_j.state.ns[4] / ndry_j;
            soil_i.state.ns[5] += n_mass * soil_j.state.ns[5] / ndry_j;
            soil_i.state.Σe += n_mass * CP_D_MOL(FT) * (soil_j.s_aux.t - T₀(FT)) / soil_i.t_aux.δz;
        end;
    end;

    # balance the air volume between top soil and atmosphere
    s_dry = top_soil.state.ns[1] + top_soil.state.ns[2] + top_soil.state.ns[4] + top_soil.state.ns[5];
    a_dry = air.state.ns[1] + air.state.ns[2] + air.state.ns[4] + air.state.ns[5];
    s_max = (air.state.p_air - saturation_vapor_pressure(top_soil.s_aux.t, top_soil.s_aux.ψ * 1000000)) * top_soil.t_aux.δz *
            max(0, top_soil.trait.vc.Θ_SAT - top_soil.state.θ - top_soil.state.θ_ice) / (GAS_R(FT) * top_soil.s_aux.t);

    # if soil air is not saturated, it can absorb more air from the atmosphere
    if s_max > s_dry
        n_mass = min(s_max - s_dry, a_dry / 2);
        top_soil.state.ns[1] += n_mass * air.state.ns[1] / a_dry;
        top_soil.state.ns[2] += n_mass * air.state.ns[2] / a_dry;
        top_soil.state.ns[4] += n_mass * air.state.ns[4] / a_dry;
        top_soil.state.ns[5] += n_mass * air.state.ns[5] / a_dry;
        top_soil.state.Σe += n_mass * CP_D_MOL(FT) * (air.s_aux.t - T₀(FT)) / top_soil.t_aux.δz;
        if !PRESCRIBE_AIR
            air.state.ns[1] -= n_mass * air.state.ns[1] / a_dry;
            air.state.ns[2] -= n_mass * air.state.ns[2] / a_dry;
            air.state.ns[4] -= n_mass * air.state.ns[4] / a_dry;
            air.state.ns[5] -= n_mass * air.state.ns[5] / a_dry;
            air.state.Σe -= n_mass * CP_D_MOL(FT) * air.s_aux.t;
        end;
    elseif s_max < s_dry
        n_mass = s_dry - s_max;
        top_soil.state.ns[1] -= n_mass * top_soil.state.ns[1] / s_dry;
        top_soil.state.ns[2] -= n_mass * top_soil.state.ns[2] / s_dry;
        top_soil.state.ns[4] -= n_mass * top_soil.state.ns[4] / s_dry;
        top_soil.state.ns[5] -= n_mass * top_soil.state.ns[5] / s_dry;
        top_soil.state.Σe -= n_mass * CP_D_MOL(FT) * (top_soil.s_aux.t - T₀(FT)) / top_soil.t_aux.δz;
        if !PRESCRIBE_AIR
            air.state.ns[1] += n_mass * top_soil.state.ns[1] / s_dry;
            air.state.ns[2] += n_mass * top_soil.state.ns[2] / s_dry;
            air.state.ns[4] += n_mass * top_soil.state.ns[4] / s_dry;
            air.state.ns[5] += n_mass * top_soil.state.ns[5] / s_dry;
            air.state.Σe += n_mass * CP_D_MOL(FT) * top_soil.s_aux.t;
        end;
    end;

    return nothing
end;
