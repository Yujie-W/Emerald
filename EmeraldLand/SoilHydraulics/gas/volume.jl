# This file contains functions to balance the air volume in the soil

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Jun-30: move function out of soil_budget!
#     2023-Jul-06: sort the order of gas volume balance and water volume balance
#     2023-Jul-06: add DEBUG code block
#     2023-Jul-06: add PRESCRIBE_AIR mode to avoid the errors due to mass balance in air
#     2023-Oct-07: limit the volume change from the source to 1/2 of the total dry air (soil and air)
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
    soils = spac.soils;

    # balance the air volume among soil layers from lower to upper layers
    air = spac.airs[1];
    for i in length(soils)-1:-1:1
        # upper layer is soil_i and lower is soil_j
        soil_i = soils[i];
        soil_j = soils[i+1];

        # compute the air moles in the lower layer
        ndry_i = soil_i.state.ns[1] + soil_i.state.ns[2] + soil_i.state.ns[4] + soil_i.state.ns[5];
        ndry_j = soil_j.state.ns[1] + soil_j.state.ns[2] + soil_j.state.ns[4] + soil_j.state.ns[5];
        nmax_j = (air.state.p_air - saturation_vapor_pressure(soil_j.auxil.t, soil_j.auxil.ψ * 1000000)) * soil_j.auxil.δz * (soil_j.state.vc.Θ_SAT - soil_j.state.θ) / (GAS_R(FT) * soil_j.auxil.t);

        # if ndry_j == nmax_j, no air needs to be transferred from/to the upper layer

        # if nmax_j > ndry_j and ndry_i > 0, air needs to be transferred from the upper layer
        if (nmax_j > ndry_j) && (ndry_i > 0)
            n_mass = min(nmax_j - ndry_j, ndry_i / 2);
            soil_j.state.ns[1] += n_mass * soil_i.state.ns[1] / ndry_i;
            soil_j.state.ns[2] += n_mass * soil_i.state.ns[2] / ndry_i;
            soil_j.state.ns[4] += n_mass * soil_i.state.ns[4] / ndry_i;
            soil_j.state.ns[5] += n_mass * soil_i.state.ns[5] / ndry_i;
            soil_j.state.Σe += n_mass * CP_D_MOL(FT) * soil_i.auxil.t / soil_j.auxil.δz;
            soil_i.state.ns[1] -= n_mass * soil_i.state.ns[1] / ndry_i;
            soil_i.state.ns[2] -= n_mass * soil_i.state.ns[2] / ndry_i;
            soil_i.state.ns[4] -= n_mass * soil_i.state.ns[4] / ndry_i;
            soil_i.state.ns[5] -= n_mass * soil_i.state.ns[5] / ndry_i;
            soil_i.state.Σe -= n_mass * CP_D_MOL(FT) * soil_i.auxil.t / soil_i.auxil.δz;
        end;

        # if nmax_j > ndry_j but upper layer is saturated, the lower layer will need to suck some water from upper layer to balance the air volume
        # compute the equilibrate mole of air from upper layer when it reaches its field capacity
        if (nmax_j > ndry_j) && (soil_i.state.θ >= soil_i.state.vc.Θ_SAT)
            θ_fc_up = soil_θ(soil_i.state.vc, -1 * ρg_MPa(FT) * soil_i.auxil.δz / relative_surface_tension(soil_i.auxil.t));
            v_mass = min((nmax_j - ndry_j) * GAS_R(FT) * soil_j.auxil.t / air.state.p_air, (soil_i.state.vc.Θ_SAT - θ_fc_up) * soil_i.auxil.δz);
            soil_j.state.θ += v_mass / soil_j.auxil.δz;
            soil_i.state.θ -= v_mass / soil_i.auxil.δz;
            soil_j.state.Σe += v_mass * ρ_H₂O(FT) * CP_L(FT) * soil_i.auxil.t / soil_j.auxil.δz;
            soil_i.state.Σe -= v_mass * ρ_H₂O(FT) * CP_L(FT) * soil_i.auxil.t / soil_i.auxil.δz;
        end;

        # if nmax_j < ndry_j and ndry_j > 0, air needs to be transferred to the upper layer (does not matter whether upper layer is saturated or not)
        if (nmax_j < ndry_j) && (ndry_j > 0)
            n_mass = ndry_j - nmax_j;
            soil_j.state.ns[1] -= n_mass * soil_j.state.ns[1] / ndry_j;
            soil_j.state.ns[2] -= n_mass * soil_j.state.ns[2] / ndry_j;
            soil_j.state.ns[4] -= n_mass * soil_j.state.ns[4] / ndry_j;
            soil_j.state.ns[5] -= n_mass * soil_j.state.ns[5] / ndry_j;
            soil_j.state.Σe -= n_mass * CP_D_MOL(FT) * soil_j.auxil.t / soil_j.auxil.δz;
            soil_i.state.ns[1] += n_mass * soil_j.state.ns[1] / ndry_j;
            soil_i.state.ns[2] += n_mass * soil_j.state.ns[2] / ndry_j;
            soil_i.state.ns[4] += n_mass * soil_j.state.ns[4] / ndry_j;
            soil_i.state.ns[5] += n_mass * soil_j.state.ns[5] / ndry_j;
            soil_i.state.Σe += n_mass * CP_D_MOL(FT) * soil_j.auxil.t / soil_i.auxil.δz;
        end;
    end;

    # balance the air volume between top soil and atmosphere
    soil = soils[1];
    s_dry = soil.state.ns[1] + soil.state.ns[2] + soil.state.ns[4] + soil.state.ns[5];
    a_dry = air.state.ns[1] + air.state.ns[2] + air.state.ns[4] + air.state.ns[5];
    s_max = (air.state.p_air - saturation_vapor_pressure(soil.auxil.t, soil.auxil.ψ * 1000000)) * soil.auxil.δz * (soil.state.vc.Θ_SAT - soil.state.θ) / (GAS_R(FT) * soil.auxil.t);

    # if soil air is not saturated, it can absorb more air from the atmosphere
    if s_max > s_dry
        n_mass = min(s_max - s_dry, a_dry / 2);
        soil.state.ns[1] += n_mass * air.state.ns[1] / a_dry;
        soil.state.ns[2] += n_mass * air.state.ns[2] / a_dry;
        soil.state.ns[4] += n_mass * air.state.ns[4] / a_dry;
        soil.state.ns[5] += n_mass * air.state.ns[5] / a_dry;
        soil.state.Σe += n_mass * CP_D_MOL(FT) * air.auxil.t / soil.auxil.δz;
        if !PRESCRIBE_AIR
            air.state.ns[1] -= n_mass * air.state.ns[1] / a_dry;
            air.state.ns[2] -= n_mass * air.state.ns[2] / a_dry;
            air.state.ns[4] -= n_mass * air.state.ns[4] / a_dry;
            air.state.ns[5] -= n_mass * air.state.ns[5] / a_dry;
            air.state.Σe -= n_mass * CP_D_MOL(FT) * air.auxil.t;
        end;
    elseif s_max < s_dry
        n_mass = s_dry - s_max;
        soil.state.ns[1] -= n_mass * soil.state.ns[1] / s_dry;
        soil.state.ns[2] -= n_mass * soil.state.ns[2] / s_dry;
        soil.state.ns[4] -= n_mass * soil.state.ns[4] / s_dry;
        soil.state.ns[5] -= n_mass * soil.state.ns[5] / s_dry;
        soil.state.Σe -= n_mass * CP_D_MOL(FT) * soil.auxil.t / soil.auxil.δz;
        if !PRESCRIBE_AIR
            air.state.ns[1] += n_mass * soil.state.ns[1] / s_dry;
            air.state.ns[2] += n_mass * soil.state.ns[2] / s_dry;
            air.state.ns[4] += n_mass * soil.state.ns[4] / s_dry;
            air.state.ns[5] += n_mass * soil.state.ns[5] / s_dry;
            air.state.Σe += n_mass * CP_D_MOL(FT) * soil.auxil.t;
        end;
    end;

    return nothing
end;
