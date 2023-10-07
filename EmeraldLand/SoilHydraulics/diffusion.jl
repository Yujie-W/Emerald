#######################################################################################################################################################################################################
#
# Changes to the function
# General
#    2023-Jun-29: copy function out of soil_budget!
#    2023-Jun-29: add code to display debug information
#     2023-Jul-06: add info into DEBUG code block
#     2023-Sep-07: add ALLOW_SOIL_EVAPORATION check
#     2023-Sep-07: fix a typo in the concentration calculations
#     2023-Sep-09: fix a typo in the concentration calculations
#
#######################################################################################################################################################################################################
"""
#
    soil_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Compute the diffusion rate among soil layers, given
- `config` SPAC configuration
- `spac` SPAC model

#
    soil_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Update diffusion rate among soil layers (and thus water and energy budgets), given
- `config` SPAC configuration
- `spac` SPAC model
- `δt` time step

"""
function soil_diffusion! end

soil_diffusion!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT} = (
    (; SOILS) = spac;

    # run the diffusion
    for soil in SOILS
        _δθ = max(0, soil.VC.Θ_SAT - soil.θ);
        if _δθ == 0
            soil.state.ns[1] = 0;
            soil.state.ns[2] = 0;
            soil.state.ns[3] = 0;
            soil.state.ns[4]  = 0;
            soil.state.ns[5]  = 0;
        else
            soil.state.ns[1] += soil.∂n∂t[1] * δt;
            soil.state.ns[2] += soil.∂n∂t[2] * δt;
            soil.state.ns[3] += soil.∂n∂t[3] * δt;
            soil.state.ns[4] += soil.∂n∂t[4] * δt;
            soil.state.ns[5] += soil.∂n∂t[5] * δt;
        end;
    end;

    return nothing
);
