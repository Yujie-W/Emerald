# This file contains functions to prescribe environmental variables and traits

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Apr-19: add the method to update the dignostic variables from air temperature
#     2022-Apr-19: add options p_CO₂, p_H₂O, rh, t, vpd, and wind (defaults are nothing)
#     2022-Nov-22: add option to change CO₂ partial pressure through ppm
#     2023-Oct-18: initialize air layer when any of the air mole fraction or temperature changes
#
#######################################################################################################################################################################################################
"""

    prescribe_air!(
                air::AirLayer{FT};
                f_CO₂::Union{Number,Nothing} = nothing,
                p_CO₂::Union{Number,Nothing} = nothing,
                p_H₂O::Union{Number,Nothing} = nothing,
                rh::Union{Number,Nothing} = nothing,
                t::Union{Number,Nothing} = nothing,
                vpd::Union{Number,Nothing} = nothing,
                wind::Union{Number,Nothing} = nothing) where {FT}

Update the environmental conditions (such as saturated vapor pressure and relative humidity) of the air surrounding the leaf, given
- `air` `AirLayer` type structure
- `f_CO₂` CO₂ concentration in `ppm`. Optional, default is nothing
- `p_CO₂` CO₂ partial pressure in `Pa`. Optional, default is nothing
- `p_H₂O` Vapor pressure in `Pa`. Optional, default is nothing
- `rh` Relatibe humidity (fraction). Optional, default is nothing
- `t` Air temperature in `K`. Optional, default is nothing
- `vpd` Vapor pressure deficit `Pa`. Optional, default is nothing
- `wind` Wind speed in `m s⁻¹`. Optional, default is nothing

"""
function prescribe_air!(
            air::AirLayer{FT};
            f_CO₂::Union{Number,Nothing} = nothing,
            p_CO₂::Union{Number,Nothing} = nothing,
            p_H₂O::Union{Number,Nothing} = nothing,
            rh::Union{Number,Nothing} = nothing,
            t::Union{Number,Nothing} = nothing,
            vpd::Union{Number,Nothing} = nothing,
            wind::Union{Number,Nothing} = nothing
) where {FT}
    if !isnothing(t) air.auxil.t = t; end;
    if !isnothing(wind) air.auxil.wind = wind; end;
    if !isnothing(f_CO₂) air.auxil.f_CO₂ = f_CO₂; air.auxil.ps[2] = air.auxil.f_CO₂ * air.state.p_air * 1e-6; end;
    if !isnothing(p_CO₂) air.auxil.ps[2] = p_CO₂; air.auxil.f_CO₂ = air.auxil.ps[2] / air.state.p_air * 1e6; end;
    if !isnothing(p_H₂O) air.auxil.ps[3] = p_H₂O; end;
    if !isnothing(rh) air.auxil.ps[3] = saturation_vapor_pressure(air.auxil.t) * rh; end;
    if !isnothing(vpd) air.auxil.ps[3] = max(0, saturation_vapor_pressure(air.auxil.t) - vpd); end;

    # if any of temperature or air mole fraction changes, re-initialize the air layer
    if !isnothing(t) || !isnothing(f_CO₂) || !isnothing(p_CO₂) || !isnothing(p_H₂O) || !isnothing(rh) || !isnothing(vpd)
        initialize_struct!(air);
    end;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-19: add function to update or prescribe swcs
#     2022-Oct-19: add option to prescribe t_soil profile
#     2023-Mar-28: update total energy in soil and leaf when prescribing swc and temperature
#     2023-Jun-12: update soil trace gas as well
#     2023-Jun-13: update N₂ and O₂ based on soil water content
#     2023-Jun-13: add soil gas energy into soil e
#     2023-Jun-15: make sure prescribed swc does not exceed the limits
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Oct-07: add 0.01 to the water vapor volume per soil layer
#     2023-Oct-18: re-initialize soil layer when any of the soil water content or temperature changes
#
#######################################################################################################################################################################################################

"""

    prescribe_soil!(spac::BulkSPAC{FT}; swcs::Union{Tuple,Nothing} = nothing, t_soils::Union{Tuple,Nothing} = nothing) where {FT}

Update the physiological parameters of the SPAC, given
- `config` Configuration for `BulkSPAC`
- `swcs` Soil water content at different layers. Optional, default is nothing
- `t_soils` Soil temperature at different layers. Optional, default is nothing

"""
function prescribe_soil!(spac::BulkSPAC{FT}; swcs::Union{Tuple,Nothing} = nothing, t_soils::Union{Tuple,Nothing} = nothing) where {FT}
    airs = spac.airs;
    soils = spac.soils;

    # prescribe soil water content (within [Θ_RES,Θ_SAT])
    if !isnothing(swcs)
        for i in eachindex(swcs)
            soils[i].state.θ = max(soils[i].state.vc.Θ_RES + eps(FT), min(soils[i].state.vc.Θ_SAT - eps(FT), swcs[i]));
            initialize_struct!(soils[i], airs[1]);
        end;
    end;

    # prescribe soil temperature
    if !isnothing(t_soils)
        for i in eachindex(t_soils)
            soils[i].auxil.t = t_soils[i];
            initialize_struct!(soils[i], airs[1]);
        end;
    end;

    # if any of soil water content or temperature changes, re-initialize the soil layer
    if !isnothing(swcs) || !isnothing(t_soils)
        for soil in soils
            initialize_struct!(soil, airs[1]);
        end;
    end;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Oct-19: add method to update or prescribe cab, car, lai, Vcmax and Jmax TD, t_leaf, vcmax profile
#     2022-Nov-21: fix a bug related to Vcmax profile (no global simulations are impacted)
#     2023-May-11: add ci to the option list
#     2023-May-19: use δlai per canopy layer
#     2023-Aug-25: add option to set up hydraulic conductance profiles for root, trunk, branches, and leaves
#     2023-Aug-27: fix a typo in the computation of k profiles (reverse the denominator and numerator)
#     2023-Oct-02: run energy initialization when LAI or t_leaf is updated
#     2023-Oct-18: recalculate canopy structural parameters when LAI, cab, car, ci is updated
#     2024-Jan-11: add option to prescribe sai
#     2024-Feb-08: fix the issue in Vcmax profiles (was opposite)
#     2024-Feb-08: add support to C3CytoState and C3VJPState
#
#######################################################################################################################################################################################################
"""

    prescribe_traits!(
                config::SPACConfiguration{FT},
                spac::BulkSPAC{FT};
                cab::Union{Number,Nothing} = nothing,
                car::Union{Number,Nothing} = nothing,
                ci::Union{Number,Nothing} = nothing,
                kmax::Union{Number,Tuple,Nothing} = nothing,
                lai::Union{Number,Nothing} = nothing,
                sai::Union{Number,Nothing} = nothing,
                swcs::Union{Tuple,Nothing} = nothing,
                t_clm::Union{Number,Nothing} = nothing,
                t_leaf::Union{Number,Nothing} = nothing,
                t_soils::Union{Tuple,Nothing} = nothing,
                vcmax::Union{Number,Nothing} = nothing,
                vcmax_expo::Union{Number,Nothing} = nothing) where {FT}

Update the physiological parameters of the SPAC, given
- `spac` Soil plant air continuum
- `config` Configuration for `BulkSPAC`
- `cab` Chlorophyll content. Optional, default is nothing
- `car` Carotenoid content. Optional, default is nothing
- `ci` Clumping index. Optional, default is nothing
- `kmax` Maximum hydraulic conductance. Optional, default is nothing
- `lai` Leaf area index. Optional, default is nothing
- `sai` Stem area index. Optional, default is nothing
- `t_clm` Moving average temperature to update Vcmax and Jmax temperature dependencies. Optional, default is nothing
- `t_leaf` Leaf temperature. Optional, default is nothing
- `vcmax` Vcmax25 at the top of canopy. Optional, default is nothing
- `vcmax_expo` Exponential tuning factor to adjust Vcmax25. Optional, default is nothing

"""
function prescribe_traits!(
            config::SPACConfiguration{FT},
            spac::BulkSPAC{FT};
            cab::Union{Number,Nothing} = nothing,
            car::Union{Number,Nothing} = nothing,
            ci::Union{Number,Nothing} = nothing,
            kmax::Union{Number,Tuple,Nothing} = nothing,
            lai::Union{Number,Nothing} = nothing,
            sai::Union{Number,Nothing} = nothing,
            t_clm::Union{Number,Nothing} = nothing,
            t_leaf::Union{Number,Nothing} = nothing,
            vcmax::Union{Number,Nothing} = nothing,
            vcmax_expo::Union{Number,Nothing} = nothing
) where {FT}
    (; T_CLM) = config;
    branches = spac.plant.branches;
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    roots = spac.plant.roots;
    sbulk = spac.soil_bulk;
    trunk = spac.plant.trunk;
    n_layer = length(leaves);

    # update chlorophyll and carotenoid contents (and spectra)
    if !isnothing(cab)
        for leaf in leaves
            leaf.bio.state.cab = cab;
        end;
    end;

    if !isnothing(car)
        for leaf in leaves
            leaf.bio.state.car = car;
        end;
    end;

    if !isnothing(cab) || !isnothing(car)
        plant_leaf_spectra!(config, spac);
    end;

    # update LAI and Vcmax (with scaling factor)
    if !isnothing(lai)
        can_str.state.lai = lai;
        can_str.state.δlai = lai .* ones(FT, n_layer) ./ n_layer;
        can_str.auxil.x_bnds = ([0; [sum(can_str.state.δlai[1:i]) + sum(can_str.state.δsai[1:i]) for i in 1:n_layer]] ./ -(lai + can_str.state.sai));
        for irt in 1:n_layer
            ilf = n_layer - irt + 1;
            leaves[ilf].xylem.state.area = sbulk.state.area * can_str.state.δlai[irt];
        end;
    end;

    if !isnothing(vcmax)
        leaves[end].photosystem.state.v_cmax25 = vcmax;
    end;

    if !isnothing(vcmax) || !isnothing(lai)
        for irt in 1:n_layer
            ilf = n_layer - irt + 1;
            ratio = isnothing(vcmax_expo) ? 1 : exp(-vcmax_expo * sum(can_str.state.δlai[1:irt-1]));
            leaf = leaves[ilf];
            if leaf.photosystem.state isa C3VJPState
                leaf.photosystem.state.v_cmax25 = leaves[end].photosystem.state.v_cmax25 * ratio;
                leaf.photosystem.state.j_max25 = leaves[end].photosystem.state.v_cmax25 * 1.67 * ratio;
                leaf.photosystem.state.r_d25 = leaves[end].photosystem.state.v_cmax25 * 0.015 * ratio;
                leaf.photosystem.auxil._t = 0;
            elseif leaf.photosystem.state isa C3CytoState
                leaf.photosystem.state.v_cmax25 = leaves[end].photosystem.state.v_cmax25 * ratio;
                leaf.photosystem.state.b₆f = leaves[end].photosystem.state.v_cmax25 * 7 / 300 * ratio;
                leaf.photosystem.state.r_d25 = leaves[end].photosystem.state.v_cmax25 * 0.015 * ratio;
                leaf.photosystem.auxil._t = 0;
            else
                error("Vcmax profile is only available for C3VJPState and C3CytoState.");
            end;
        end;
    end;

    # update CI
    if !isnothing(ci)
        can_str.state.Ω_A = ci;
        can_str.state.Ω_B = 0;
        can_str.auxil.ci = ci;
    end;

    # update SAI
    if !isnothing(sai)
        can_str.state.sai = sai;
        can_str.state.δsai = sai .* ones(FT, n_layer) ./ n_layer;
        can_str.auxil.x_bnds = ([0; [sum(can_str.state.δlai[1:i]) + sum(can_str.state.δsai[1:i]) for i in 1:n_layer]] ./ -(can_str.state.lai + sai));
    end;

    # update Vcmax and Jmax TD
    if !isnothing(t_clm)
        for leaf in leaves
            if T_CLM
                leaf.photosystem.state.TD_VCMAX.ΔSV = 668.39 - 1.07 * (t_clm - T₀(FT));
                leaf.photosystem.state.TD_JMAX.ΔSV = 659.70 - 0.75 * (t_clm - T₀(FT));
            end;
        end;
    end;

    # update kmax
    if !isnothing(kmax)
        # set up the kmax assuming 50% resistance in root, 25% in stem, and 25% in leaves
        ks = if kmax isa Number
            trunk_percent = trunk.xylem.state.Δh / (trunk.xylem.state.Δh + branches[end].xylem.state.Δh);
            (2 * kmax, 4 * kmax / trunk_percent, 4 * kmax / (1 - trunk_percent), 4 * kmax)
        else
            @assert length(kmax) == 4 "kmax must be a number or a tuple of length 4";
            kmax
        end;

        # partition kmax into the roots based on xylem area
        for root in roots
            # root.xylem.state.k_max = root.xylem.state.area / trunk.xylem.state.area * ks[1] * root.xylem.state.l / root.xylem.state.area;
            root.xylem.state.k_max = ks[1] * root.xylem.state.l / trunk.xylem.state.area;
        end;
        trunk.xylem.state.k_max = ks[2] * trunk.xylem.state.l / trunk.xylem.state.area;
        for stem in branches
            #stem.xylem.state.kmax = stem.xylem.state.area / trunk.xylem.state.area * ks[3] * stem.xylem.state.l / stem.xylem.state.area;
            stem.xylem.state.kmax = ks[3] * stem.xylem.state.l / trunk.xylem.state.area;
        end;
        for leaf in leaves
            leaf.xylem.state.k_max = ks[4] / (can_str.state.lai * sbulk.state.area);
        end;
    end;

    # prescribe leaf temperature
    if !isnothing(t_leaf)
        for leaf in leaves
            leaf.energy.auxil.t = t_leaf;
        end;
    end;

    # re-initialize leaf energy if LAI or t_leaf is updated
    if !isnothing(lai) || !isnothing(t_leaf)
        for leaf in leaves
            initialize_struct!(leaf; k_sla = leaf.xylem.state.k_max);
        end;
    end;

    # recalibrate the canopy structure parameters if any of LAI, Cab, Car, etc. is updated
    if !isnothing(lai) || !isnothing(cab) || !isnothing(car) || !isnothing(ci)
        canopy_structure!(config, spac);
    end;

    return nothing
end;
