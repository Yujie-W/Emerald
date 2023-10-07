# This function contains the soil budgets for water, gas, and energy

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-07: add function to run the soil water condensation or evaporation
#
#######################################################################################################################################################################################################
"""

    soil_water_condensation!(soil::SoilLayer{FT}) where {FT}

Run the soil water condensation or evaporation, given
- `soil` `SoilLayer` type soil layer

"""
function soil_water_condensation!(soil::SoilLayer{FT}) where {FT}
    p_sat = saturation_vapor_pressure(soil.auxil.t, soil.auxil.ψ * 1000000);
    n_con = soil.state.ns[3] - p_sat * max(0, soil.state.vc.Θ_SAT - soil.state.θ) * soil.auxil.δz / (GAS_R(FT) * soil.auxil.t);
    v_liq = n_con * M_H₂O(FT) / ρ_H₂O(FT);

    soil.auxil.n_con = n_con;
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

    soil_budgets!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Run the soil budgets for water and gas, given
- `config` `SPACConfiguration` type configuration
- `spac` `MultiLayerSPAC` type SPAC
- `δt` time step

"""
function soil_budgets!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}
    (; SOIL_BULK, SOILS) = spac;

    for i in eachindex(SOILS)
        # run the soil gas budget
        for j in 1:5
            SOILS[i].state.ns[j] += SOILS[i].auxil.∂n∂t[j] * δt;
        end;

        # run the water transport (mass flow)
        SOILS[i].state.θ += SOILS[i].auxil.∂θ∂t * δt;

        # run the soil water condensation or evaporation
        soil_water_condensation!(SOILS[i]);
    end;

    volume_balance!(config, spac);

    # compute surface runoff
    if SOILS[1].state.θ > SOILS[1].state.vc.Θ_SAT
        SOIL_BULK.auxil.runoff = (SOILS[1].state.θ - SOILS[1].state.vc.Θ_SAT) * SOILS[1].auxil.δz * ρ_H₂O(FT) / M_H₂O(FT);;
        SOILS[1].state.θ = SOILS[1].state.vc.Θ_SAT;
    end;

    return nothing
end;