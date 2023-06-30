#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-13: add function for water budget
#     2022-Jun-14: use K_MAX and ΔZ and remove K_REF
#     2022-Jun-14: rescale rain for layer 1
#     2022-Jun-14: use METEO.rain
#     2022-Jun-14: add function for soil energy budget
#     2022-Jun-14: use METEO.rain and METEO.t_precip
#     2022-Jun-14: add net radiation energy to top soil
#     2022-Jun-15: add controller to make sure soil layers do not over saturate
#     2022-Jun-15: merge the soil_water! and soil_energy! to soil_budget!
#     2022-Jun-16: move time stepper controller to SoilPlantAirContinuum.jl
#     2022-Jul-26: fix the unit of rain, mass flow, and root extraction (all in mol s⁻¹)
#     2022-Sep-07: allow soil water oversaturation
#     2023-Mar-27: fix a typo when updating e per layer (should use ΔZ per layer rather than the first layer)
#     2023-Apr-07: fix a typo when updating water content in saturated soil layers
#     2023-Apr-08: make runoff a cumulative value within a time interval
#     2023-Jun-13: add trace gas diffusions
#     2023-Jun-13: add diffusion related water and energy budgets
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#
#######################################################################################################################################################################################################
"""

    soil_budget!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}) where {FT<:AbstractFloat}

Update the marginal increase of soil water content and energy per layer, given
- `spac` `MultiLayerSPAC` SPAC
- `config` Configuration for `MultiLayerSPAC`

"""
soil_budget!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}) where {FT<:AbstractFloat} = (
    soil_mass_flow!(config, spac);
    soil_diffusion!(config, spac);
    soil_source_sink!(spac);

    return nothing
);


"""

    soil_mass_flow!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Compute the water and energy budget related to mass flow, given
- `config` Configuration for `MultiLayerSPAC`
- `spac` `MultiLayerSPAC` SPAC

"""
function soil_mass_flow!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}
    (; METEO, SOIL) = spac;
    (; DIM_SOIL) = config;
    LAYERS = SOIL.LAYERS;

    # update soil k, kd, ψ, and λ_thermal for each soil layer (0.5 for tortuosity factor)
    for _slayer in LAYERS
        _slayer.k          = relative_hydraulic_conductance(_slayer.VC, _slayer.θ) * _slayer.VC.K_MAX * relative_viscosity(_slayer.t) / _slayer.ΔZ;
        _slayer.ψ          = soil_ψ_25(_slayer.VC, _slayer.θ; oversaturation = true) * relative_surface_tension(_slayer.t);
        _slayer._kd        = 0.5 * max(0, _slayer.VC.Θ_SAT - _slayer.θ) / _slayer.ΔZ;
        _slayer._λ_thermal = (_slayer.Λ_THERMAL + _slayer.θ * Λ_THERMAL_H₂O(FT)) / _slayer.ΔZ;
        _slayer.∂e∂t       = 0;
        _slayer.∂n∂t      .= 0;
        _slayer.∂θ∂t       = 0;
    end;

    # update k, δψ, and flow rate among layers
    LAYERS[1].∂θ∂t += METEO.rain * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[1].ΔZ;
    LAYERS[1].∂e∂t += METEO.rain * CP_L_MOL(FT) * METEO.t_precip;
    LAYERS[1].∂e∂t += SOIL.ALBEDO.r_net_lw + SOIL.ALBEDO.r_net_sw;
    for _i in 1:DIM_SOIL-1
        SOIL._k[_i]         = 1 / (2 / LAYERS[_i].k + 2 / LAYERS[_i+1].k);
        SOIL._δψ[_i]        = LAYERS[_i].ψ - LAYERS[_i+1].ψ + ρg_MPa(FT) * (LAYERS[_i].Z - LAYERS[_i+1].Z);
        SOIL._q[_i]         = SOIL._k[_i] * SOIL._δψ[_i];
        SOIL._λ_thermal[_i] = 1 / (2 / LAYERS[_i]._λ_thermal + 2 / LAYERS[_i+1]._λ_thermal);
        SOIL._δt[_i]        = LAYERS[_i].t - LAYERS[_i+1].t;
        SOIL._q_thermal[_i] = SOIL._λ_thermal[_i] * SOIL._δt[_i];

        # if flow into the lower > 0, but the lower layer is already saturated, set the flow to 0
        if (SOIL._q[_i] > 0) && (LAYERS[_i+1].θ >= LAYERS[_i+1].VC.Θ_SAT)
            SOIL._q[_i] = 0;
        end;

        # if flow into the lower < 0, but the upper layer is already saturated, set the flow to 0
        if (SOIL._q[_i] < 0) && (LAYERS[_i].θ >= LAYERS[_i].VC.Θ_SAT)
            SOIL._q[_i] = 0;
        end;

        # if both layers are oversaturated, move the oversaturated part from lower layer to upper layer
        if (LAYERS[_i].θ >= LAYERS[_i].VC.Θ_SAT) && (LAYERS[_i+1].θ > LAYERS[_i+1].VC.Θ_SAT)
            SOIL._q[_i] = -1 * (LAYERS[_i+1].θ - LAYERS[_i+1].VC.Θ_SAT) * LAYERS[_i+1].ΔZ * ρ_H₂O(FT) / M_H₂O(FT);
        end;

        LAYERS[_i  ].∂θ∂t -= SOIL._q[_i] * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[_i].ΔZ;
        LAYERS[_i+1].∂θ∂t += SOIL._q[_i] * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[_i+1].ΔZ;
        LAYERS[_i  ].∂e∂t -= SOIL._q_thermal[_i];
        LAYERS[_i+1].∂e∂t += SOIL._q_thermal[_i];
        LAYERS[_i  ].∂e∂t -= SOIL._q[_i] * CP_L_MOL(FT) * LAYERS[_i].t;
        LAYERS[_i+1].∂e∂t += SOIL._q[_i] * CP_L_MOL(FT) * LAYERS[_i].t;
    end;

    return nothing
end





"""

    soil_source_sink!(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}

Update the source/sink terms for the soil layers, given
- `spac` the SPAC model

"""
function soil_source_sink!(spac::MultiLayerSPAC{FT}) where {FT<:AbstractFloat}
    (; ROOTS, ROOTS_INDEX, SOIL) = spac;
    LAYERS = SOIL.LAYERS;

    # loop through the roots and compute the source/sink terms
    for _i in eachindex(ROOTS)
        LAYERS[ROOTS_INDEX[_i]].∂θ∂t -= root_sink(ROOTS[_i]) * M_H₂O(FT) / ρ_H₂O(FT) / SOIL.AREA / LAYERS[ROOTS_INDEX[_i]].ΔZ;
        LAYERS[ROOTS_INDEX[_i]].∂e∂t -= root_sink(ROOTS[_i]) / SOIL.AREA * CP_L_MOL(FT) * LAYERS[_i].t;
    end;

    return nothing
end
