#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Apr-19: separate this function as an individual step of the SPAC module (1st step)
#     2022-Jul-12: rename function to update!
#
#######################################################################################################################################################################################################
"""
This function updates the environmental conditions and the soil-plant-air-continuum. Supported functionalities are for
- AirLayer
- SoilPlantAirContinuum

"""
function update! end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Apr-19: add the method to update the dignostic variables from air temperature
#     2022-Apr-19: add options p_CO₂, p_H₂O, rh, t, vpd, and wind (defaults are nothing)
#     2022-Apr-19: update docs and history log
#     2022-Jul-12: rename function to update!
#     2022-Jul-12: remove FT control to options
#     2022-Oct-19: add method to update or prescribe cab, car, lai, swcs, Vcmax and Jmax TD, t_leaf, vcmax profile
#     2022-Oct-19: air.rh and air.p_H₂O_sat have been removed in an earlier version ClimaCache
#     2022-Oct-19: add method to prescribe t_soil profile
#     2022-Nov-21: fix a bug related to Vcmax profile (no global simulations are impacted)
#     2022-Nov-22: add option to change CO₂ partial pressure through ppm
#     2023-Mar-28: fix a typo when updating t_soil
#     2023-Mar-28: update total energy in soil and leaf when prescribing swc and temperature
#     2023-May-11: add ci to the option list
#     2023-May-19: use δlai per canopy layer
#     2023-Jun-12: update soil trace gas as well
#     2023-Jun-12: fix trace gas initialization
#     2023-Jun-13: update N₂ and O₂ based on soil water content
#     2023-Jun-13: add soil gas energy into soil e
#     2023-Jun-15: make sure prescribed swc does not exceed the limits
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Aug-25: add option to set yo hydraulic conductance profiles for root, trunk, branches, and leaves
#     2023-Aug-27: fix a typo in the computation of k profiles (reverse the denominator and numerator)
#     2023-Sep-07: add ALLOW_LEAF_CONDENSATION and T_CLM checks
#     2023-Oct-02: run energy initialization when LAI or t_leaf is updated
#     2023-Oct-07: add 0.01 to the water vapor volume per soil layer
#
#######################################################################################################################################################################################################

"""

    update!(config::SPACConfiguration{FT},
            spac::BulkSPAC{FT},;
            cab::Union{Number,Nothing} = nothing,
            car::Union{Number,Nothing} = nothing,
            ci::Union{Number,Nothing} = nothing,
            kmax::Union{Number,Tuple,Nothing} = nothing,
            lai::Union{Number,Nothing} = nothing,
            swcs::Union{Tuple,Nothing} = nothing,
            t_clm::Union{Number,Nothing} = nothing,
            t_leaf::Union{Number,Nothing} = nothing,
            t_soils::Union{Tuple,Nothing} = nothing,
            vcmax::Union{Number,Nothing} = nothing,
            vcmax_expo::Union{Number,Nothing} = nothing
    ) where {FT}

Update the physiological parameters of the SPAC, given
- `spac` Soil plant air continuum
- `config` Configuration for `BulkSPAC`
- `cab` Chlorophyll content. Optional, default is nothing
- `car` Carotenoid content. Optional, default is nothing
- `ci` Clumping index. Optional, default is nothing
- `kmax` Maximum hydraulic conductance. Optional, default is nothing
- `lai` Leaf area index. Optional, default is nothing
- `swcs` Soil water content at different layers. Optional, default is nothing
- `t_clm` Moving average temperature to update Vcmax and Jmax temperature dependencies. Optional, default is nothing
- `t_leaf` Leaf temperature. Optional, default is nothing
- `t_soils` Soil temperature at different layers. Optional, default is nothing
- `vcmax` Vcmax25 at the top of canopy. Optional, default is nothing
- `vcmax_expo` Exponential tuning factor to adjust Vcmax25. Optional, default is nothing

"""
update!(config::SPACConfiguration{FT},
        spac::BulkSPAC{FT};
        cab::Union{Number,Nothing} = nothing,
        car::Union{Number,Nothing} = nothing,
        ci::Union{Number,Nothing} = nothing,
        kmax::Union{Number,Tuple,Nothing} = nothing,
        lai::Union{Number,Nothing} = nothing,
        swcs::Union{Tuple,Nothing} = nothing,
        t_clm::Union{Number,Nothing} = nothing,
        t_leaf::Union{Number,Nothing} = nothing,
        t_soils::Union{Tuple,Nothing} = nothing,
        vcmax::Union{Number,Nothing} = nothing,
        vcmax_expo::Union{Number,Nothing} = nothing
) where {FT} = (
    (; DIM_LAYER, T_CLM) = config;
    airs = spac.airs;
    branches = spac.plant.branches;
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    roots = spac.plant.roots;
    sbulk = spac.soil_bulk;
    soils = spac.soils;
    trunk = spac.plant.trunk;

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
        can_str.state.δlai = lai .* ones(FT, DIM_LAYER) ./ DIM_LAYER;
        can_str.auxil.x_bnds = ([0; [sum(can_str.state.δlai[1:i]) + sum(can_str.state.δsai[1:i]) for i in 1:DIM_LAYER]] ./ -(lai + can_str.state.sai));
        for i in 1:DIM_LAYER
            leaves[i].xylem.state.area = sbulk.state.area * can_str.state.δlai[i];
        end;

        # make sure leaf area index setup and energy are correct
        for i in eachindex(leaves)
            leaves[i].xylem.state.area = sbulk.state.area * can_str.state.δlai[i];
            initialize_struct!(leaves[i]);
        end;
    end;
    if !isnothing(vcmax)
        leaves[1].photosystem.state.v_cmax25 = vcmax;
    end;
    if !isnothing(vcmax) || !isnothing(lai)
        for i in 2:DIM_LAYER
            _scaling = isnothing(vcmax_expo) ? 1 : exp(-vcmax_expo * sum(can_str.state.δlai[1:i-1]));
            leaves[i].photosystem.state.v_cmax25 = leaves[1].photosystem.state.v_cmax25 * _scaling;
            leaves[i].photosystem.state.j_max25 = leaves[1].photosystem.state.v_cmax25 * 1.67 * _scaling;
            leaves[i].photosystem.state.r_d25 = leaves[1].photosystem.state.v_cmax25 * 0.015 * _scaling;
            leaves[i].photosystem.auxil._t = 0;
        end;
    end;

    # update CI
    if !isnothing(ci)
        can_str.auxil.ci = ci;
        can_str.state.Ω_A = ci;
        can_str.state.Ω_B = 0;
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
        _ks = if kmax isa Number
            _trunk_percent = trunk.HS.ΔH / (trunk.HS.ΔH + branches[end].HS.ΔH);
            (2 * kmax, 4 * kmax / _trunk_percent, 4 * kmax / (1 - _trunk_percent), 4 * kmax)
        else
            @assert length(kmax) == 4 "kmax must be a number or a tuple of length 4";
            kmax
        end;

        # partition kmax into the roots based on xylem area
        for root in roots
            root.HS.K_X = root.HS.AREA / trunk.HS.AREA * _ks[1] * root.HS.L / root.HS.AREA;
        end;
        trunk.HS.K_X = _ks[2] * trunk.HS.L / trunk.HS.AREA;
        for stem in branches
            stem.HS.K_X = stem.HS.AREA / trunk.HS.AREA * _ks[3] * stem.HS.L / stem.HS.AREA;
        end;
        for leaf in leaves
            leaf.HS.K_SLA = _ks[4] / (can_str.state.lai * sbulk.state.area);
        end;
    end;

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

    # prescribe leaf temperature
    if !isnothing(t_leaf)
        for leaf in leaves
            leaf.energy.auxil.t = t_leaf;
            leaf.energy.auxil.cp = heat_capacitance(leaf);
            leaf.energy.state.Σe = leaf.energy.auxil.cp * leaf.energy.auxil.t;
        end;

        # make sure leaf area index setup and energy are correct
        for i in eachindex(leaves)
            leaves[i].xylem.state.area = sbulk.state.area * can_str.state.δlai[i];
            initialize_struct!(leaves[i]);
        end;
    end;

    return nothing
);
