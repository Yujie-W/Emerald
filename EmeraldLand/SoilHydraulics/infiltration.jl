#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Jun-29: copy function out of soil_budget!
#     2023-Jun-30: add DEBUG mode to check for NaNs
#     2023-Jul-06: add info into DEBUG code block
#     2023-Sep-07: add ALLOW_SOIL_EVAPORATION check
#     2023-Sep-07: add integrators for soil water budget
#
#######################################################################################################################################################################################################
"""

    soil_infiltration!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Update soil water content and energy per layer, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC
- `δt` Time step

"""
function soil_infiltration! end

soil_infiltration!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT} = (
    (; SOILS) = spac;

    # run the water transport (condensation + mass flow)
    for i in eachindex(SOILS)
        soil = SOILS[i];

        # account for evaporation and condensation to/from the air space
        _ps = saturation_vapor_pressure(soil.t, soil.auxil.ψ * 1000000);
        _δθ_v = (soil.state.ns[3] / soil.auxil.δz - _ps * max(0, soil.state.vc.Θ_SAT - soil.θ) / (GAS_R(FT) * soil.t)) * M_H₂O(FT) / ρ_H₂O(FT);

        soil.state.θ += _δθ_v;
        soil.state.Σe += _δθ_v * ρ_H₂O(FT) * CP_V(FT) * soil.t; # this energy is transferred from/to air, so use CP_V
        soil.state.Σe += _δθ_v * ρ_H₂O(FT) * latent_heat_vapor(soil.t);

        # account for mass flow
        soil.state.θ += soil.auxil.∂θ∂t * δt;
        soil.state.Σe += soil.auxil.∂e∂t * δt / soil.auxil.δz;
    end;

    return nothing
);
