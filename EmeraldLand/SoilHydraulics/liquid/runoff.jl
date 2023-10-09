# This function contains the function to compute the runoff of the soil

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-09: add function to compute the runoff of the soil from the lowest layer to the top layer
#
#######################################################################################################################################################################################################
"""

    soil_water_runoff!(spac::MultiLayerSPAC{FT}) where {FT}

Compute the runoff of the soil from the lowest layer to the top layer, given
- `spac` `MultiLayerSPAC` type SPAC

"""
function soil_water_runoff!(spac::MultiLayerSPAC{FT}) where {FT}
    (; SOIL_BULK, SOILS) = spac;

    # iterate from the lowest soil layer to the top
    N = length(SOILS);
    for i in N-1:-1:1
        soili = SOILS[i];
        soilj = SOILS[i+1];

        # if the lower layer is oversaturated, the oversaturated water will flow to the upper layer (as well as energy)
        if soilj.state.θ > soilj.state.vc.Θ_SAT
            # compute the water flow rate and energy associated with the flow
            v_flow = (soilj.state.θ - soilj.state.vc.Θ_SAT) * soilj.auxil.δz;
            e_flow = v_flow * ρ_H₂O(FT) * CP_L(FT) * soilj.auxil.t;

            # update the water and energy in the layers
            soili.state.θ  += v_flow / soili.auxil.δz;
            soilj.state.θ  -= v_flow / soilj.auxil.δz;
            soili.state.Σe += e_flow / soili.auxil.δz;
            soilj.state.Σe -= e_flow / soilj.auxil.δz;
        end;
    end;

    # run the soil water runoffn from the top layer (energy will be updated in EnergyBudgets.jl because of the need to call heat_capacitance function)
    if SOILS[1].state.θ > SOILS[1].state.vc.Θ_SAT
        SOIL_BULK.auxil.runoff = (SOILS[1].state.θ - SOILS[1].state.vc.Θ_SAT) * SOILS[1].auxil.δz * ρ_H₂O(FT) / M_H₂O(FT);
        SOILS[1].state.θ = SOILS[1].state.vc.Θ_SAT;
    end;

    return nothing
end;