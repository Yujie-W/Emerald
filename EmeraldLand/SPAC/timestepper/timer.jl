# This file contains function to calculate the time step size

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-15: add function to make sure the soil layers do not over saturate or drain
#     2022-Jun-18: add controller for soil and leaf temperatures
#     2022-Aug-18: add option θ_on to enable/disable soil water budget
#     2022-Aug-31: add controller for leaf stomatal conductance
#     2022-Sep-07: remove soil oversaturation controller, and add a Δθ <= 0.01 controller
#     2022-Oct-22: add option t_on to enable/disable soil and leaf energy budgets
#     2023-Mar-27: add controller for trunk and branch temperatures
#     2023-Sep-30: add time controller of junction water
#     2023-Oct-18: add time controller of junction temperature
#     2024-Feb-29: add methods for soil and plant time controllers (for cleaner reading experience, may need to add air in the future)
#     2024-Jul-24: make sure the water content of the junction does not exceed the minimum water content
#     2024-Aug-06: add junction pressure change controller (not change more than 0.1 MPa per time step)
#     2024-Sep-03: run stem and branck T check as well (was deactivated by accident)
#
#######################################################################################################################################################################################################
"""

        adjusted_time(spac::BulkSPAC{FT}, δt::FT) where {FT}

Return adjusted time that soil does not over saturate or drain, given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC
- `δt` Time step

"""
function adjusted_time end;

adjusted_time(spac::BulkSPAC{FT}, δt::FT) where {FT} = (
    new_δt::FT = δt;

    # adjust the time based on soil
    for soil in spac.soils
        new_δt = adjusted_time(soil, new_δt);
    end;

    # adjust the time based on plant
    new_δt = adjusted_time(spac.plant, spac.canopy.structure.trait.lai, new_δt);

    return new_δt
);

adjusted_time(soil::SoilLayer{FT}, δt::FT) where {FT} = (
    new_δt::FT = δt;

    # make sure soil water content will not change more than 0.01 m³ m⁻³ per time step
    new_δt = min(FT(0.01) / abs(soil.auxil.∂θ∂t), new_δt);
    if isnan(new_δt) || new_δt < 0.001 <= δt
        @error "NaN or very small δt detected when adjusting δt based on soil moisture" soil.auxil.∂θ∂t;
        return error("NaN detected in adjusted_time!")
    end;

    # make sure soil water will not completely drain (to Θ_RES)
    if soil.auxil.∂θ∂t < 0
        new_δt = min((soil.trait.vc.Θ_RES - soil.state.θ) / soil.auxil.∂θ∂t, new_δt);
        if isnan(new_δt) || new_δt < 0.001 <= δt
            @error "NaN or very small δt detected that could drain the soil" soil.trait.vc.Θ_RES soil.state.θ soil.auxil.∂θ∂t;
            return error("NaN detected in adjusted_time")
        end;
    end;

    # make sure soil temperature does not change more than 1 K per time step
    ∂T∂t = soil.auxil.∂e∂t / soil.s_aux.cp;
    new_δt = min(1 / abs(∂T∂t), new_δt);
    if isnan(new_δt) || new_δt < 0.001 <= δt
        @error "NaN or very small δt detected when adjusting δt based on soil temperature" soil.auxil.∂e∂t soil.s_aux.cp ∂T∂t;
        return error("NaN detected in adjusted_time")
    end;

    return new_δt
);

adjusted_time(plant::Plant{FT}, lai::FT, δt::FT) where {FT} = adjusted_time(plant, lai, δt, plant.leaves[1]);

adjusted_time(plant::Plant{FT}, lai::FT, δt::FT, ::CanopyLayer{FT}) where {FT} = (
    new_δt::FT = δt;

    # make sure root temperature does not change more than 1 K per time step
    for root in plant.roots
        ∂T∂t = root.energy.auxil.∂e∂t / root.energy.s_aux.cp;
        new_δt = min(1 / abs(∂T∂t), new_δt);
        if isnan(new_δt) || new_δt < 0.001 <= δt
            @error "NaN or very small δt detected when adjusting δt based on root temperature" root.energy.auxil.∂e∂t root.energy.s_aux.cp ∂T∂t;
            return error("NaN detected in adjusted_time")
        end;
    end;

    # make sure junction temperature does not change more than 1 K per time step
    ∂T∂t = plant.junction.auxil.∂e∂t / plant.junction.s_aux.cp;
    new_δt = min(1 / abs(∂T∂t), new_δt);
    if isnan(new_δt) || new_δt < 0.001 <= δt
        @error "NaN or very small δt detected when adjusting δt based on junction temperature" plant.junction.auxil.∂e∂t plant.junction.s_aux.cp ∂T∂t;
        return error("NaN detected in adjusted_time")
    end;

    # make sure junction water does not change more than 10 mol per time step
    new_δt = min(10 / abs(plant.junction.auxil.∂w∂t), new_δt);
    if isnan(new_δt) || new_δt < 0.001 <= δt
        @error "NaN or very small δt detected when adjusting δt based on junction water" plant.junction.auxil.∂w∂t;
        return error("NaN detected in adjusted_time")
    end;

    # make sure the water content of the junction does not exceed the minimum water content (half through here)
    new_δt = min((plant.junction.state.v_storage - plant.junction.trait.v_max * plant.junction.trait.pv.residual) / abs(plant.junction.auxil.∂w∂t) / 2, new_δt);
    if isnan(new_δt) || new_δt < 0.001 <= δt
        @error "NaN or very small δt detected when adjusting δt based on junction water storage" plant.junction.state.v_storage plant.junction.auxil.∂w∂t;
        return error("NaN detected in adjusted_time")
    end;

    # make sure that the junction pressure does not change more than 0.1 MPa per time step (change the order to avoid new_v < 0)
    new_v = plant.junction.state.v_storage + plant.junction.auxil.∂w∂t * new_δt;
    new_p = capacitance_pressure(plant.junction.trait.pv, new_v / plant.junction.trait.v_max, plant.junction.s_aux.t);
    if abs(new_p - plant.junction.s_aux.pressure) > 0.1
        new_δt *= FT(0.1) / abs(new_p - plant.junction.s_aux.pressure);
    end;
    if isnan(new_δt) || new_δt < 0.001 <= δt
        @error "NaN or very small δt detected when adjusting δt based on junction pressure change" plant.junction.s_aux.pressure new_p;
        return error("NaN detected in adjusted_time")
    end;

    # make sure trunk temperature does not change more than 1 K per time step
    ∂T∂t = plant.trunk.energy.auxil.∂e∂t / plant.trunk.energy.s_aux.cp;
    new_δt = min(1 / abs(∂T∂t), new_δt);
    if isnan(new_δt) || new_δt < 0.001 <= δt
        @error "NaN or very small δt detected when adjusting δt based on trunk temperature" plant.trunk.energy.auxil.∂e∂t plant.trunk.energy.s_aux.cp ∂T∂t;
        return error("NaN detected in adjusted_time")
    end;

    # make sure each branch temperature does not change more than 1 K per time step
    for stem in plant.branches
        ∂T∂t = stem.energy.auxil.∂e∂t / stem.energy.s_aux.cp;
        new_δt = min(1 / abs(∂T∂t), new_δt);
        if isnan(new_δt) || new_δt < 0.001 <= δt
            @error "NaN or very small δt detected when adjusting δt based on branch temperature" stem.energy.auxil.∂e∂t stem.energy.s_aux.cp ∂T∂t;
            return error("NaN detected in adjusted_time")
        end;
    end;

    # leaves adjustments are required only when LAI > 0
    if lai <= 0
        return new_δt
    end;

    # make sure each leaf temperature does not change more than 1 K per time step
    for leaf in plant.leaves
        ∂T∂t = leaf.energy.auxil.∂e∂t / leaf.energy.s_aux.cp;
        new_δt = min(1 / abs(∂T∂t), new_δt);
        if isnan(new_δt) || new_δt < 0.001 <= δt
            @error "NaN or very small δt detected when adjusting δt based on leaf temperature" leaf.energy.auxil.∂e∂t leaf.energy.s_aux.cp ∂T∂t;
            return error("NaN detected in adjusted_time")
        end;
    end;

    # make sure each leaf stomatal conductances do not change more than 0.01 mol m⁻² s⁻¹
    for leaf in plant.leaves
        for ∂g∂t in leaf.flux.auxil.∂g∂t
            new_δt = min(FT(0.01) / abs(∂g∂t), new_δt);
            if isnan(new_δt) || new_δt < 0.001 <= δt
                @error "NaN or very small δt detected when adjusting δt based on leaf stomatal conductance dYdt" ∂g∂t;
                return error("NaN detected in adjusted_time")
            end;
        end;
    end;

    return new_δt
);
