module SoilHydraulics

using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak

using ..Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, M_H₂O, Λ_THERMAL_H₂O, ρ_H₂O, ρg_MPa
using ..Namespace: MultiLayerSPAC, NonSteadyStateFlow, Root, SPACConfiguration, SteadyStateFlow, VanGenuchten
using ..PhysicalChemistry: diffusive_coefficient, latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure

import ..Namespace: BrooksCorey


#######################################################################################################################################################################################################
#
# Changes made to this constructor
# General
#     2021-Sep-30: move this function out of BrooksCorey struct as an external method for the constructor (avoid dependency on ConstrainedRootSolvers)
#     2022-Apr-19: fix documentation
#     2022-Jul-15: BrooksCorey field changed, modify the constructor accordingly
#
#######################################################################################################################################################################################################
"""

    BrooksCorey{FT}(vg::VanGenuchten{FT}) where {FT<:AbstractFloat}

A constructor for BrooksCorey to create BrooksCorey type soil from VanGenuchten type, given
- `vg` `VanGenuchten` type soil water retention curve
"""
BrooksCorey{FT}(vg::VanGenuchten{FT}) where {FT<:AbstractFloat} = (
    _bc = BrooksCorey{FT}(K_MAX = vg.K_MAX, B = 1, TYPE = vg.TYPE, Ψ_SAT = 0.001, Θ_SAT = vg.Θ_SAT, Θ_RES = vg.Θ_RES);

    # generate data to fit
    _Θs   = range(vg.Θ_RES+FT(1e-2); stop=vg.Θ_SAT-FT(1e-2), length=30);
    _Ψ_vG = -1 .* soil_ψ_25.([vg], _Θs);
    _Ψ_BC = similar(_Ψ_vG);
    _Ψ_DF = similar(_Ψ_vG);

    # function to fit BrooksCorey parameters
    @inline _fit(x) = (
        _bc.B = x[1];
        _bc.Ψ_SAT = x[2];
        _Ψ_BC .= -1 .* soil_ψ_25.([_bc], _Θs);
        _Ψ_DF .= (log.(_Ψ_BC) .- log.(_Ψ_vG)) .^ 2;
        return -sum(_Ψ_DF)
    );

    _st  = SolutionToleranceND{FT}([1e-3, 1e-6], 30);
    _ms  = ReduceStepMethodND{FT}(x_mins = FT[1e-3, 1e-6], x_maxs = FT[ 100, 1000], x_inis = [(2*vg.N-1) / (vg.N-1), 1 / (vg.α)], Δ_inis = FT[0.1, 1e-3]);
    _sol = find_peak(_fit, _ms, _st);
    _bc.B = _sol[1];
    _bc.Ψ_SAT = _sol[2];

    return _bc
);


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Sep-30: create this function to work with two soil types using either VanGenuchten or BrooksCorey function
#     2022-Sep-07: add option to allow for soil water oversaturation
#     2022-Oct-20: set a minimum effective θ at eps(FT)
#
#######################################################################################################################################################################################################
"""

    soil_ψ_25(bc::BrooksCorey{FT}, θ::FT; oversaturation::Bool = false) where {FT<:AbstractFloat}
    soil_ψ_25(vg::VanGenuchten{FT}, θ::FT; oversaturation::Bool = false) where {FT<:AbstractFloat}

Return the soil metric potential, given
- `bc` or `vg` `BrooksCorey` or `VanGenuchten` type structure
- `θ` Soil volumetric water content (absolute value)
- `oversaturation` If true, allow for soil water oversaturation

"""
function soil_ψ_25 end

soil_ψ_25(bc::BrooksCorey{FT}, θ::FT; oversaturation::Bool = false) where {FT<:AbstractFloat} = (
    (; B, Ψ_SAT, Θ_RES, Θ_SAT) = bc;

    # calculate effective θ
    _θ_e = max(eps(FT), (θ - Θ_RES) / (Θ_SAT - Θ_RES));

    # if _θ_e >= 1, use segmented function
    if _θ_e >= 1
        return oversaturation ? (θ - Θ_RES) - Ψ_SAT : -Ψ_SAT
    end;

    return -Ψ_SAT / (_θ_e ^ B)
);

soil_ψ_25(vg::VanGenuchten{FT}, θ::FT; oversaturation::Bool = false) where {FT<:AbstractFloat} = (
    (; M, N, α, Θ_RES, Θ_SAT) = vg;

    # calculate effective θ
    _θ_e = max(eps(FT), (θ - Θ_RES) / (Θ_SAT - Θ_RES));

    # if _θ_e >= 1, use segmented function
    if _θ_e >= 1
        return oversaturation ? θ - Θ_RES : FT(0)
    end;

    return -1 * (_θ_e ^ (-1/M) - 1) ^ (1/N) / α
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-01: migrate function to version v0.3 of PlantHydraulics v0.2
#
#######################################################################################################################################################################################################
"""

    soil_θ(bc::BrooksCorey{FT}, ψ_25::FT) where {FT<:AbstractFloat}
    soil_θ(vg::VanGenuchten{FT}, ψ_25::FT) where {FT<:AbstractFloat}

Return the soil water content, given
- `bc` or `vg` `BrooksCorey` or `VanGenuchten` type structure
- `ψ_25` Soil metric potential corrected to 25 Celcius

"""
function soil_θ end

soil_θ(bc::BrooksCorey{FT}, ψ_25::FT) where {FT<:AbstractFloat} = (
    (; B, Ψ_SAT, Θ_RES, Θ_SAT) = bc;

    if ψ_25 >= 0
        return Θ_SAT
    end;

    return (-Ψ_SAT/ψ_25) ^ (1/B) * (Θ_SAT - Θ_RES) + Θ_RES
);

soil_θ(vg::VanGenuchten{FT}, ψ_25::FT) where {FT<:AbstractFloat} = (
    (; M, N, α, Θ_RES, Θ_SAT) = vg;

    if ψ_25 >= 0
        return Θ_SAT
    end;

    return ( 1 / ( 1 + (-ψ_25 * α) ^ N ) ) ^ M * (Θ_SAT - Θ_RES) + Θ_RES
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Feb-01: migrate function to version v0.3 of PlantHydraulics
#     2022-Feb-01: add documentation
#     2022-Feb-01: rename the function to relative_hydraulic_conductance
#     2022-Feb-01: remove viscosity correction
#     2022-Apr-19: move function to SoilHydraulics from PlantHydraulics (will be imported in PlantHydraulics)
#     2022-May-31: add a controller to ψ_25 to make avoid numerical issue
#     2022-Oct-20: set a minimum relative k at eps(FT)
#
#######################################################################################################################################################################################################
"""

    relative_hydraulic_conductance(bc::BrooksCorey{FT}, θ::FT) where {FT<:AbstractFloat}
    relative_hydraulic_conductance(bc::BrooksCorey{FT}, ψ::Bool, ψ_25::FT) where {FT<:AbstractFloat}
    relative_hydraulic_conductance(vg::VanGenuchten{FT}, θ::FT) where {FT<:AbstractFloat}
    relative_hydraulic_conductance(vg::VanGenuchten{FT}, ψ::Bool, ψ_25::FT) where {FT<:AbstractFloat}

Return the relative hydraulic conductance of the soil, given
- `bc` or `vg` `BrooksCorey` or `VanGenuchten` type structure
- `θ` Soil volumetric water content (absolute value)
- `ψ` Bool to indicate that next parameter is potential
- `ψ_25` Soil metric potential at a reference temperature of 25 °C

"""
function relative_hydraulic_conductance end

relative_hydraulic_conductance(bc::BrooksCorey{FT}, θ::FT) where {FT<:AbstractFloat} = (
    (; B, Θ_RES, Θ_SAT) = bc;

    _θ_e = min(1, max(eps(FT), (θ - Θ_RES) / (Θ_SAT - Θ_RES)));

    return max(eps(FT), _θ_e ^ (2 * B + 3))
);

relative_hydraulic_conductance(bc::BrooksCorey{FT}, ψ::Bool, ψ_25::FT) where {FT<:AbstractFloat} = (
    (; B, Ψ_SAT) = bc;

    # if the potential > 0, return 1
    if ψ_25 >= 0
        return FT(1)
    end;

    _θ_e = (-Ψ_SAT / ψ_25) ^ (1 / B);

    return max(eps(FT), _θ_e ^ (2 * B + 3))
);

relative_hydraulic_conductance(vg::VanGenuchten{FT}, θ::FT) where {FT<:AbstractFloat} = (
    (; M, N, Θ_RES, Θ_SAT) = vg;

    _θ_e = min(1, max(eps(FT), (θ - Θ_RES) / (Θ_SAT - Θ_RES)));

    return max(eps(FT), sqrt(_θ_e) * (1 - (1 - _θ_e ^ (1 / M)) ^ M)^2)
);

relative_hydraulic_conductance(vg::VanGenuchten{FT}, ψ::Bool, ψ_25::FT) where {FT<:AbstractFloat} = (
    (; M, N, α) = vg;

    # if the potential > 0, return 1
    if ψ_25 >= 0
        return FT(1)
    end;

    _θ_e = (1 / (1 + (-ψ_25 * α) ^ N)) ^ M;

    return max(eps(FT), sqrt(_θ_e) * (1 - (1 - _θ_e ^ (1 / M)) ^ M)^2)
);


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
This function have two major features:
- Compute the marginal change of soil water content and energy
- Update soil water content and energy without over-saturating or draining the soil

"""
function soil_budget! end


"""

    soil_budget!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}) where {FT<:AbstractFloat}

Update the marginal increase of soil water content and energy per layer, given
- `spac` `MultiLayerSPAC` SPAC
- `config` Configuration for `MultiLayerSPAC`

"""
soil_budget!(spac::MultiLayerSPAC{FT}, config::SPACConfiguration{FT}) where {FT<:AbstractFloat} = (
    (; AIR, METEO, ROOTS, ROOTS_INDEX, SOIL) = spac;
    (; DIM_SOIL, TRACE_AIR, TRACE_CH₄, TRACE_CO₂, TRACE_H₂O, TRACE_N₂, TRACE_O₂) = config;
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

    # update diffusion rate among layers
    _alayer = AIR[1];
    _factor = 2 * LAYERS[1]._kd;
    _v_gas = max(0, LAYERS[1].VC.Θ_SAT - LAYERS[1].θ);
    LAYERS[1].∂n∂t[1] -= _factor * diffusive_coefficient(LAYERS[1].t, TRACE_CH₄, TRACE_AIR) * (LAYERS[1].TRACES.n_CH₄ / (LAYERS[1].ΔZ * _v_gas) - _alayer.p_CH₄ / (GAS_R(FT) * _alayer.t));
    LAYERS[1].∂n∂t[2] -= _factor * diffusive_coefficient(LAYERS[1].t, TRACE_CO₂, TRACE_AIR) * (LAYERS[1].TRACES.n_CO₂ / (LAYERS[1].ΔZ * _v_gas) - _alayer.p_CO₂ / (GAS_R(FT) * _alayer.t));
    LAYERS[1].∂n∂t[3] -= _factor * diffusive_coefficient(LAYERS[1].t, TRACE_H₂O, TRACE_AIR) * (LAYERS[1].TRACES.n_H₂O / (LAYERS[1].ΔZ * _v_gas) - _alayer.p_H₂O / (GAS_R(FT) * _alayer.t));
    LAYERS[1].∂n∂t[4] -= _factor * diffusive_coefficient(LAYERS[1].t, TRACE_N₂ , TRACE_AIR) * (LAYERS[1].TRACES.n_N₂  / (LAYERS[1].ΔZ * _v_gas) - _alayer.p_N₂  / (GAS_R(FT) * _alayer.t));
    LAYERS[1].∂n∂t[5] -= _factor * diffusive_coefficient(LAYERS[1].t, TRACE_O₂ , TRACE_AIR) * (LAYERS[1].TRACES.n_O₂  / (LAYERS[1].ΔZ * _v_gas) - _alayer.p_O₂  / (GAS_R(FT) * _alayer.t));

    for _i in 1:DIM_SOIL-1
        # gas diffusion
        _ratei1 = diffusive_coefficient(LAYERS[_i  ].t, TRACE_CH₄, TRACE_AIR);
        _ratei2 = diffusive_coefficient(LAYERS[_i  ].t, TRACE_CO₂, TRACE_AIR);
        _ratei3 = diffusive_coefficient(LAYERS[_i  ].t, TRACE_H₂O, TRACE_AIR);
        _ratei4 = diffusive_coefficient(LAYERS[_i  ].t, TRACE_N₂ , TRACE_AIR);
        _ratei5 = diffusive_coefficient(LAYERS[_i  ].t, TRACE_O₂ , TRACE_AIR);
        _ratej1 = diffusive_coefficient(LAYERS[_i+1].t, TRACE_CH₄, TRACE_AIR);
        _ratej2 = diffusive_coefficient(LAYERS[_i+1].t, TRACE_CO₂, TRACE_AIR);
        _ratej3 = diffusive_coefficient(LAYERS[_i+1].t, TRACE_H₂O, TRACE_AIR);
        _ratej4 = diffusive_coefficient(LAYERS[_i+1].t, TRACE_N₂ , TRACE_AIR);
        _ratej5 = diffusive_coefficient(LAYERS[_i+1].t, TRACE_O₂ , TRACE_AIR);
        _ratio1 = 2 * LAYERS[_i]._kd * _ratei1 * LAYERS[_i+1]._kd * _ratej1 / (LAYERS[_i]._kd * _ratei1 + LAYERS[_i+1]._kd * _ratej1);
        _ratio2 = 2 * LAYERS[_i]._kd * _ratei2 * LAYERS[_i+1]._kd * _ratej2 / (LAYERS[_i]._kd * _ratei2 + LAYERS[_i+1]._kd * _ratej2);
        _ratio3 = 2 * LAYERS[_i]._kd * _ratei3 * LAYERS[_i+1]._kd * _ratej3 / (LAYERS[_i]._kd * _ratei3 + LAYERS[_i+1]._kd * _ratej3);
        _ratio4 = 2 * LAYERS[_i]._kd * _ratei4 * LAYERS[_i+1]._kd * _ratej4 / (LAYERS[_i]._kd * _ratei4 + LAYERS[_i+1]._kd * _ratej4);
        _ratio5 = 2 * LAYERS[_i]._kd * _ratei5 * LAYERS[_i+1]._kd * _ratej5 / (LAYERS[_i]._kd * _ratei5 + LAYERS[_i+1]._kd * _ratej5);
        _drate1 = _ratio1 * (LAYERS[_i].TRACES.n_CH₄ / LAYERS[_i].ΔZ - LAYERS[_i+1].TRACES.n_CH₄ / LAYERS[_i+1].ΔZ);
        _drate2 = _ratio2 * (LAYERS[_i].TRACES.n_CO₂ / LAYERS[_i].ΔZ - LAYERS[_i+1].TRACES.n_CO₂ / LAYERS[_i+1].ΔZ);
        _drate3 = _ratio3 * (LAYERS[_i].TRACES.n_H₂O / LAYERS[_i].ΔZ - LAYERS[_i+1].TRACES.n_H₂O / LAYERS[_i+1].ΔZ);
        _drate4 = _ratio4 * (LAYERS[_i].TRACES.n_N₂  / LAYERS[_i].ΔZ - LAYERS[_i+1].TRACES.n_N₂  / LAYERS[_i+1].ΔZ);
        _drate5 = _ratio5 * (LAYERS[_i].TRACES.n_O₂  / LAYERS[_i].ΔZ - LAYERS[_i+1].TRACES.n_O₂  / LAYERS[_i+1].ΔZ);
        LAYERS[_i  ].∂n∂t[1] -= _drate1;
        LAYERS[_i  ].∂n∂t[2] -= _drate2;
        LAYERS[_i  ].∂n∂t[3] -= _drate3;
        LAYERS[_i  ].∂n∂t[4] -= _drate4;
        LAYERS[_i  ].∂n∂t[5] -= _drate5;
        LAYERS[_i+1].∂n∂t[1] += _drate1;
        LAYERS[_i+1].∂n∂t[2] += _drate2;
        LAYERS[_i+1].∂n∂t[3] += _drate3;
        LAYERS[_i+1].∂n∂t[4] += _drate4;
        LAYERS[_i+1].∂n∂t[5] += _drate5;

        # energy transfer related to gas diffusion
        _δe_gas = ((_drate1 + _drate2 + _drate4 + _drate5) * CP_D_MOL(FT) + _drate3 * CP_V_MOL(FT)) * LAYERS[_i].t;
        LAYERS[_i  ].∂e∂t -= _δe_gas;
        LAYERS[_i+1].∂e∂t += _δe_gas;
    end;

    # loop through the roots and compute the source/sink terms
    for _i in eachindex(ROOTS)
        LAYERS[ROOTS_INDEX[_i]].∂θ∂t -= root_sink(ROOTS[_i]) * M_H₂O(FT) / ρ_H₂O(FT) / SOIL.AREA / LAYERS[ROOTS_INDEX[_i]].ΔZ;
        LAYERS[ROOTS_INDEX[_i]].∂e∂t -= root_sink(ROOTS[_i]) / SOIL.AREA * CP_L_MOL(FT) * LAYERS[_i].t;
    end;

    return nothing
);


"""

    soil_budget!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT<:AbstractFloat}

Run soil water and energy budget, given
- `spac` `MultiLayerSPAC` SPAC
- `δt` Time step

"""
soil_budget!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT<:AbstractFloat} = (
    (; AIR, SOIL) = spac;
    LAYERS = SOIL.LAYERS;

    # run the diffusion
    for _slayer in LAYERS
        _δv = _slayer.ΔZ * max(0, _slayer.VC.Θ_SAT - _slayer.θ);
        if _δv == 0
            _slayer.TRACES.n_CH₄ = 0;
            _slayer.TRACES.n_CO₂ = 0;
            _slayer.TRACES.n_H₂O = 0;
            _slayer.TRACES.n_N₂  = 0;
            _slayer.TRACES.n_O₂  = 0;
        else
            _slayer.TRACES.n_CH₄ += _slayer.∂n∂t[1] * δt;
            _slayer.TRACES.n_CO₂ += _slayer.∂n∂t[2] * δt;
            _slayer.TRACES.n_H₂O += _slayer.∂n∂t[3] * δt;
            _slayer.TRACES.n_N₂  += _slayer.∂n∂t[4] * δt;
            _slayer.TRACES.n_O₂  += _slayer.∂n∂t[5] * δt;
        end;
    end;

    # compute air volume change using ideal gas law (total energy change accordingly)
    _alayer = AIR[1];
    for _i in length(LAYERS):-1:1
        _slayer = LAYERS[_i];
        _n_air = (_alayer.P_AIR - saturation_vapor_pressure(_slayer.t, _slayer.ψ * 1000000)) * _slayer.ΔZ / (GAS_R(FT) * _slayer.t);
        _n_dry = _slayer.TRACES.n_CH₄ + _slayer.TRACES.n_CO₂ + _slayer.TRACES.n_N₂ + _slayer.TRACES.n_O₂;
        _n_rat = (_n_air - _n_dry) / _n_dry;
        _slayer.TRACES.n_CH₄ += _n_rat * _slayer.TRACES.n_CH₄;
        _slayer.TRACES.n_CO₂ += _n_rat * _slayer.TRACES.n_CO₂;
        _slayer.TRACES.n_N₂  += _n_rat * _slayer.TRACES.n_N₂;
        _slayer.TRACES.n_O₂  += _n_rat * _slayer.TRACES.n_O₂;
        _slayer.e += (_n_air - _n_dry) * CP_D_MOL(FT) * _slayer.t;
        if _i == 1
            _alayer.n_CH₄ -= _n_rat * _slayer.TRACES.n_CH₄;
            _alayer.n_CO₂ -= _n_rat * _slayer.TRACES.n_CO₂;
            _alayer.n_N₂  -= _n_rat * _slayer.TRACES.n_N₂;
            _alayer.n_O₂  -= _n_rat * _slayer.TRACES.n_O₂;
            _alayer.e -= (_n_air - _n_dry) * CP_D_MOL(FT) * _slayer.t;
        else
            _tlayer = LAYERS[_i-1];
            _tlayer.TRACES.n_CH₄ -= _n_rat * _slayer.TRACES.n_CH₄;
            _tlayer.TRACES.n_CO₂ -= _n_rat * _slayer.TRACES.n_CO₂;
            _tlayer.TRACES.n_N₂  -= _n_rat * _slayer.TRACES.n_N₂;
            _tlayer.TRACES.n_O₂  -= _n_rat * _slayer.TRACES.n_O₂;
            _tlayer.e -= (_n_air - _n_dry) * CP_D_MOL(FT) * _slayer.t;
        end;
    end;

    # run the water transport (condensation + mass flow)
    for _slayer in LAYERS
        # account for evaporation and condensation to/from the air space
        _ps = saturation_vapor_pressure(_slayer.t, _slayer.ψ * 1000000);
        _δθ_v = (_slayer.TRACES.n_H₂O / _slayer.ΔZ - _ps * max(0, _slayer.VC.Θ_SAT - _slayer.θ) / (GAS_R(FT) * _slayer.t)) * M_H₂O(FT) / ρ_H₂O(FT);
        _slayer.θ += _δθ_v;
        _slayer.e += _δθ_v * ρ_H₂O(FT) * CP_L(FT) * _slayer.t;
        _slayer.e += _δθ_v * ρ_H₂O(FT) * latent_heat_vapor(_slayer.t);

        # account for mass flow
        _slayer.θ += _slayer.∂θ∂t * δt;
        _slayer.e += _slayer.∂e∂t * δt / _slayer.ΔZ;
    end;

    # compute surface runoff
    if LAYERS[1].θ > LAYERS[1].VC.Θ_SAT
        # compute top soil temperature and top soil energy out due to runoff
        _cp = LAYERS[1].CP * LAYERS[1].ρ + LAYERS[1].θ * ρ_H₂O(FT) * CP_L(FT);
        _t  = LAYERS[1].e / _cp;
        _runoff = (LAYERS[1].θ - LAYERS[1].VC.Θ_SAT) * LAYERS[1].ΔZ * ρ_H₂O(FT) / M_H₂O(FT);
        LAYERS[1].θ = LAYERS[1].VC.Θ_SAT;
        LAYERS[1].e -= _runoff / LAYERS[1].ΔZ * CP_L_MOL(FT) * _t;
        SOIL.runoff += _runoff;
    end;

    # update soil temperature at each layer (top layer t will be same as _t above)
    for _slayer in LAYERS
        _cp_gas = (_slayer.TRACES.n_H₂O * CP_V_MOL(FT) + (_slayer.TRACES.n_CH₄ + _slayer.TRACES.n_CO₂ + _slayer.TRACES.n_N₂ + _slayer.TRACES.n_O₂) * CP_D_MOL(FT)) / _slayer.ΔZ;
        _slayer._cp = _slayer.ρ * _slayer.CP + _slayer.θ * ρ_H₂O(FT) * CP_L(FT) + _cp_gas;
        _slayer.t = _slayer.e / _slayer._cp;
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-13: add a utility function to read root water sink
#
#######################################################################################################################################################################################################
"""

    root_sink(root::Root{FT}) where {FT<:AbstractFloat}

Return root water update, given
- `root` `Root` type struct that may contain non- and steady state flow

"""
function root_sink end

root_sink(root::Root{FT}) where {FT<:AbstractFloat} = root_sink(root.HS.FLOW);

root_sink(mode::SteadyStateFlow{FT}) where {FT<:AbstractFloat} = mode.flow;

root_sink(mode::NonSteadyStateFlow{FT}) where {FT<:AbstractFloat} = mode.f_in;


end # module
