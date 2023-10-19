# This function contains the soil budgets for water, gas, and energy

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-07: add function to run the soil water condensation or evaporation
#     2023-Oct-07: add 0.01 to the water vapor volume per soil layer
#     2023-Oct-09: remove condensation or evaporation of water vapor from the soil trace gas moles
#
#######################################################################################################################################################################################################
"""

    soil_water_condensation!(soil::SoilLayer{FT}) where {FT}

Run the soil water condensation or evaporation, given
- `soil` `SoilLayer` type soil layer

"""
function soil_water_condensation!(soil::SoilLayer{FT}) where {FT}
    p_sat = saturation_vapor_pressure(soil.auxil.t, soil.auxil.ψ * 1000000);
    n_con = soil.state.ns[3] - p_sat * (max(0, soil.state.vc.Θ_SAT - soil.state.θ) * soil.auxil.δz + FT(0.01)) / (GAS_R(FT) * soil.auxil.t);
    v_liq = n_con * M_H₂O(FT) / ρ_H₂O(FT);

    soil.auxil.n_con = n_con;
    soil.state.ns[3] -= n_con;
    soil.state.θ += v_liq / soil.auxil.δz;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-07: add function to run the soil budgets for water and gas
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

        # run the soil water condensation or evaporation
        soil_water_condensation!(soil);
    end;

    # compute surface runoff and volume balance for the air
    soil_water_runoff!(spac);
    volume_balance!(config, spac);

    return nothing
end;
