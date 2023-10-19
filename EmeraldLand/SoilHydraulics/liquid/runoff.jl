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
    N = length(soils);
    for i in N-1:-1:1
        soili = soils[i];
        soilj = soils[i+1];

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
    top_soil = soils[1];
    if top_soil.state.θ > top_soil.state.vc.Θ_SAT
        sbulk.auxil.runoff = (top_soil.state.θ - top_soil.state.vc.Θ_SAT) * top_soil.auxil.δz * ρ_H₂O(FT) / M_H₂O(FT);
        top_soil.state.θ = top_soil.state.vc.Θ_SAT;
    end;

    return nothing
end;
