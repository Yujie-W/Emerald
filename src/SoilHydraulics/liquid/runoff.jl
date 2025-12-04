# This function contains the function to compute the runoff of the soil

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-09: add function to compute the runoff of the soil from the lowest layer to the top layer
#
#######################################################################################################################################################################################################
"""

    soil_water_runoff!(spac::BulkSPAC{FT}) where {FT}

Compute the runoff of the soil from the lowest layer to the top layer, given
- `spac` `BulkSPAC` type SPAC

"""
function soil_water_runoff!(spac::BulkSPAC{FT}) where {FT}
    sbulk = spac.soil_bulk;
    soils = spac.soils;

    # iterate from the lowest soil layer to the top
    # TODO: ice volume is not accounted for, so the water flow may not be correct
    N = length(soils);
    for i in N-1:-1:1
        soili = soils[i];
        soilj = soils[i+1];

        # if the lower layer is oversaturated, the oversaturated water will flow to the upper layer (as well as energy)
        if soilj.state.θ > soilj.trait.vc.Θ_SAT
            # compute the water flow rate and energy associated with the flow
            v_flow = (soilj.state.θ - soilj.trait.vc.Θ_SAT) * soilj.t_aux.δz;
            e_flow = v_flow * ρ_H₂O(FT) * CP_L(FT) * soilj.s_aux.t;

            # update the water and energy in the layers
            soili.state.θ  += v_flow / soili.t_aux.δz;
            soilj.state.θ  -= v_flow / soilj.t_aux.δz;
            soili.state.Σe += e_flow / soili.t_aux.δz;
            soilj.state.Σe -= e_flow / soilj.t_aux.δz;
        end;
    end;

    # run the soil water runoff from the top layer (energy will be updated in EnergyBudgets.jl because of the need to call heat_capacitance function)
    # TODO: ice volume is not accounted for, so the water flow may not be correct
    top_soil = soils[1];
    if top_soil.state.θ > top_soil.trait.vc.Θ_SAT
        sbulk.auxil.runoff = (top_soil.state.θ - top_soil.trait.vc.Θ_SAT) * top_soil.t_aux.δz * ρ_H₂O(FT) / M_H₂O(FT);
        top_soil.state.θ = top_soil.trait.vc.Θ_SAT;
    end;

    return nothing
end;
