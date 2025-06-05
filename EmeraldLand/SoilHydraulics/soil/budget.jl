# This function contains the soil budgets for water, gas, and energy

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-07: add function to run the soil water condensation or evaporation
#     2023-Oct-07: add 0.01 to the water vapor volume per soil layer
#     2023-Oct-09: remove condensation or evaporation of water vapor from the soil trace gas moles
#     2025-Jun-05: account for ice volume in the calculation of condensation or evaporation
#
#######################################################################################################################################################################################################
"""

    soil_water_condensation!(soil::SoilLayer{FT}) where {FT}

Run the soil water condensation or evaporation, given
- `soil` `SoilLayer` type soil layer

"""
function soil_water_condensation!(soil::SoilLayer{FT}) where {FT}
    p_sat = saturation_vapor_pressure(soil.s_aux.t, soil.s_aux.ψ * 1000000);
    n_con = soil.state.ns[3] - p_sat * (max(0, soil.trait.vc.Θ_SAT - soil.state.θ - soil.state.θ_ice) * soil.t_aux.δz + FT(0.01)) / (GAS_R(FT) * soil.s_aux.t);
    v_liq = n_con * M_H₂O(FT) / ρ_H₂O(FT);

    soil.auxil.n_con = n_con;
    soil.state.ns[3] -= n_con;
    soil.state.θ += v_liq / soil.t_aux.δz;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2025-Jun-05: add function to run the soil water freeze/thaw
#
#######################################################################################################################################################################################################
"""

    soil_water_freeze_thaw!(soil::SoilLayer{FT}) where {FT}

Run the soil water freeze or thaw, given
- `soil` `SoilLayer` type soil layer

"""
function soil_water_freeze_thaw!(soil::SoilLayer{FT}) where {FT}
    # when there is water and ice, water thaw occurs all the time when Σe > 0
    #     if the energy is enough to melt the ice, then melt all
    #     otherwise, melt as much as possible
    if soil.state.θ_ice > 0 && soil.state.Σe > 0
        e_melt = soil.state.θ_ice * ρ_H₂O(FT) * latent_heat_melt(soil.s_aux.t);
        if soil.state.Σe >= e_melt
            soil.state.θ += soil.state.θ_ice;
            soil.state.θ_ice = 0;
            soil.state.Σe -= e_melt;
        else
            dθ_ice = soil.state.Σe / e_melt * soil.state.θ_ice;
            soil.state.θ += dθ_ice;
            soil.state.θ_ice -= dθ_ice;
            soil.state.Σe -= dθ_ice / soil.state.θ_ice * e_melt;
        end;

        return nothing
    end;

    # when there is water only, make sure the energy is enough to lower the total temperature to -2 Celcius
    # then add a 0.001 fraction of ice core, so that the part below can be triggered
    if soil.state.θ_ice == 0 && soil.state.Σe < -2 * soil.s_aux.cp
        dθ_ice = min(FT(0.001), soil.state.θ - soil.trait.vc.Θ_RES);
        soil.state.θ -= dθ_ice;
        soil.state.θ_ice += dθ_ice;
        soil.state.Σe += dθ_ice * ρ_H₂O(FT) * latent_heat_melt(soil.s_aux.t);
    end;

    # when there is water and ice, water freeze occurs when Σe < 0
    #     if the energy is enough to freeze the water, then freeze all water except the residual water
    #     otherwise, freeze as much as possible
    if soil.state.θ > soil.trait.vc.Θ_RES && soil.state.θ_ice > 0 && soil.state.Σe < 0
        e_freeze = (soil.state.θ - soil.trait.vc.Θ_RES) * ρ_H₂O(FT) * latent_heat_melt(soil.s_aux.t);
        if soil.state.Σe <= -e_freeze
            dθ_ice = soil.state.θ - soil.trait.vc.Θ_RES;
            soil.state.θ -= dθ_ice;
            soil.state.θ_ice += dθ_ice;
            soil.state.Σe += e_freeze;
        else
            dθ_ice = -soil.state.Σe / e_freeze * (soil.state.θ - soil.trait.vc.Θ_RES);
            soil.state.θ -= dθ_ice;
            soil.state.θ_ice += dθ_ice;
            soil.state.Σe += dθ_ice / (soil.state.θ - soil.trait.vc.Θ_RES) * e_freeze;
        end;

        return nothing
    end;

    # do nothing otherwise
    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-07: add function to run the soil budgets for water and gas
#     2025-Jun-05: run soil freezing/thawing control at the end of each sub-step
#
#######################################################################################################################################################################################################
"""

    soil_budgets!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, δt::FT) where {FT}

Run the soil budgets for water and gas, given
- `config` `SPACConfiguration` type configuration
- `spac` `BulkSPAC` type SPAC
- `δt` time step

"""
function soil_budgets!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, δt::FT) where {FT}
    soils = spac.soils;

    # per soil layer, run the gas diffusion, mass flow, and condensation or evaporation budgets
    for soil in soils
        # run the soil gas budget
        for j in 1:5
            (soil).state.ns[j] += (soil).auxil.∂n∂t[j] * δt;
        end;

        # run the water transport (mass flow)
        (soil).state.θ += (soil).auxil.∂θ∂t * δt;

        # run the soil water condensation/evaporation and freeze/thaw
        soil_water_condensation!(soil);
        soil_water_freeze_thaw!(soil);
    end;

    # compute surface runoff and volume balance for the air
    soil_water_runoff!(spac);
    volume_balance!(config, spac);

    return nothing
end;
