module SoilHydraulics

using ..EmeraldMath.Solver: ReduceStepMethodND, SolutionToleranceND, find_peak

using ..Constant: CP_L, CP_L_MOL, M_H₂O, Λ_THERMAL_H₂O, ρ_H₂O, ρg_MPa
using ..Namespace: MonoMLTreeSPAC, NonSteadyStateFlow, Root, SteadyStateFlow, VanGenuchten
using ..PhysicalChemistry: relative_surface_tension, relative_viscosity

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
#
#######################################################################################################################################################################################################
"""
This function have two major features:
- Compute the marginal change of soil water content and energy
- Update soil water content and energy without over-saturating or draining the soil

"""
function soil_budget! end


"""

    soil_budget!(spac::MonoMLTreeSPAC{FT}) where {FT<:AbstractFloat}

Update the marginal increase of soil water content and energy per layer, given
- `spac` `MonoMLTreeSPAC` SPAC

"""
soil_budget!(spac::MonoMLTreeSPAC{FT}) where {FT<:AbstractFloat} = (
    (; METEO, ROOTS, ROOTS_INDEX, SOIL) = spac;
    LAYERS = SOIL.LAYERS;

    # update soil k, ψ, and λ_thermal for each soil layer
    for _i in 1:SOIL.DIM_SOIL
        LAYERS[_i].k          = relative_hydraulic_conductance(LAYERS[_i].VC, LAYERS[_i].θ) * LAYERS[_i].VC.K_MAX * relative_viscosity(LAYERS[_i].t) / LAYERS[_i].ΔZ;
        LAYERS[_i].ψ          = soil_ψ_25(LAYERS[_i].VC, LAYERS[_i].θ; oversaturation = true) * relative_surface_tension(LAYERS[_i].t);
        LAYERS[_i]._λ_thermal = (LAYERS[_i].Λ_THERMAL + LAYERS[_i].θ * Λ_THERMAL_H₂O(FT)) / LAYERS[_i].ΔZ;
        LAYERS[_i].∂e∂t       = 0;
        LAYERS[_i].∂θ∂t       = 0;
    end;

    # update k, δψ, and flow rate among layers
    LAYERS[1].∂θ∂t += METEO.rain * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[1].ΔZ;
    LAYERS[1].∂e∂t += METEO.rain * CP_L_MOL(FT) * METEO.t_precip;
    LAYERS[1].∂e∂t += SOIL.ALBEDO.r_net_lw + SOIL.ALBEDO.r_net_sw;
    for _i in 1:SOIL.DIM_SOIL-1
        SOIL._k[_i]         = 1 / (2 / LAYERS[_i].k + 2 / LAYERS[_i+1].k);
        SOIL._δψ[_i]        = LAYERS[_i].ψ - LAYERS[_i+1].ψ + ρg_MPa(FT) * (LAYERS[_i].Z - LAYERS[_i+1].Z);
        SOIL._q[_i]         = SOIL._k[_i] * SOIL._δψ[_i];
        SOIL._λ_thermal[_i] = 1 / (2 / LAYERS[_i]._λ_thermal + 2 / LAYERS[_i+1]._λ_thermal);
        SOIL._δt[_i]        = LAYERS[_i].t - LAYERS[_i+1].t;
        SOIL._q_thermal[_i] = SOIL._λ_thermal[_i] * SOIL._δt[_i];

        # if flow into the lower > 0, but the lower layer is already saturated, set the flow to 0
        if (SOIL._q[_i] > 0) && (SOIL.LAYERS[_i+1].θ >= SOIL.LAYERS[_i+1].VC.Θ_SAT)
            SOIL._q[_i] = 0;
        end;

        # if flow into the lower < 0, but the upper layer is already saturated, set the flow to 0
        if (SOIL._q[_i] < 0) && (SOIL.LAYERS[_i].θ >= SOIL.LAYERS[_i].VC.Θ_SAT)
            SOIL._q[_i] = 0;
        end;

        # if both layers are oversaturated, move the oversaturated part from lower layer to upper layer
        if (SOIL.LAYERS[_i].θ >= SOIL.LAYERS[_i].VC.Θ_SAT) && (SOIL.LAYERS[_i+1].θ > SOIL.LAYERS[_i+1].VC.Θ_SAT)
            SOIL._q[_i] = (SOIL.LAYERS[_i+1].θ - SOIL.LAYERS[_i+1].VC.Θ_SAT) * LAYERS[_i+1].ΔZ * ρ_H₂O(FT) / M_H₂O(FT);
        end;

        LAYERS[_i  ].∂θ∂t -= SOIL._q[_i] * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[_i].ΔZ;
        LAYERS[_i+1].∂θ∂t += SOIL._q[_i] * M_H₂O(FT) / ρ_H₂O(FT) / LAYERS[_i+1].ΔZ;
        LAYERS[_i  ].∂e∂t -= SOIL._q_thermal[_i];
        LAYERS[_i+1].∂e∂t += SOIL._q_thermal[_i];
        LAYERS[_i  ].∂e∂t -= SOIL._q[_i] * CP_L_MOL(FT) * LAYERS[_i].t;
        LAYERS[_i+1].∂e∂t += SOIL._q[_i] * CP_L_MOL(FT) * LAYERS[_i].t;
    end;

    # loop through the roots and compute the source/sink terms
    for _i in eachindex(ROOTS)
        LAYERS[ROOTS_INDEX[_i]].∂θ∂t -= root_sink(ROOTS[_i]) * M_H₂O(FT) / ρ_H₂O(FT) / SOIL.AREA / LAYERS[ROOTS_INDEX[_i]].ΔZ;
        LAYERS[ROOTS_INDEX[_i]].∂e∂t -= root_sink(ROOTS[_i]) / SOIL.AREA * CP_L_MOL(FT) * LAYERS[_i].t;
    end;

    return nothing
);


"""

    soil_budget!(spac::MonoMLTreeSPAC{FT}, δt::FT) where {FT<:AbstractFloat}

Run soil water and energy budget, given
- `spac` `MonoMLTreeSPAC` SPAC
- `δt` Time step

"""
soil_budget!(spac::MonoMLTreeSPAC{FT}, δt::FT) where {FT<:AbstractFloat} = (
    (; SOIL) = spac;

    # run the time step
    for _i in 1:SOIL.DIM_SOIL
        SOIL.LAYERS[_i].θ += SOIL.LAYERS[_i].∂θ∂t * δt;
        SOIL.LAYERS[_i].e += SOIL.LAYERS[_i].∂e∂t * δt / SOIL.LAYERS[1].ΔZ;
    end;

    # compute surface runoff
    SOIL.runoff = 0;
    if SOIL.LAYERS[1].θ > SOIL.LAYERS[1].VC.Θ_SAT
        # compute top soil temperature and top soil energy out due to runoff
        _cp = SOIL.LAYERS[1].CP * SOIL.LAYERS[1].ρ + SOIL.LAYERS[1].θ * ρ_H₂O(FT) * CP_L(FT) + SOIL.runoff / SOIL.LAYERS[1].ΔZ * CP_L_MOL(FT);
        _t  = SOIL.LAYERS[1].e / _cp;
        SOIL.runoff = (SOIL.LAYERS[1].θ - SOIL.LAYERS[1].VC.Θ_SAT) * SOIL.LAYERS[1].ΔZ * ρ_H₂O(FT) / M_H₂O(FT);
        SOIL.LAYERS[1].θ = SOIL.LAYERS[1].VC.Θ_SAT;
        SOIL.LAYERS[1].e -= SOIL.runoff / SOIL.LAYERS[1].ΔZ * CP_L_MOL(FT) * _t;
    end;

    # update soil temperature at each layer (top layer t will be same as _t above)
    for _i in 1:SOIL.DIM_SOIL
        SOIL.LAYERS[_i]._cp = SOIL.LAYERS[_i].CP * SOIL.LAYERS[_i].ρ + SOIL.LAYERS[_i].θ * ρ_H₂O(FT) * CP_L(FT);
        SOIL.LAYERS[_i].t  = SOIL.LAYERS[_i].e / SOIL.LAYERS[_i]._cp;
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
