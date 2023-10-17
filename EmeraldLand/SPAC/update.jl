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

    update!(air::AirLayer{FT};
            p_CO₂::Union{Number,Nothing} = nothing,
            p_H₂O::Union{Number,Nothing} = nothing,
            rh::Union{Number,Nothing} = nothing,
            t::Union{Number,Nothing} = nothing,
            vpd::Union{Number,Nothing} = nothing,
            wind::Union{Number,Nothing} = nothing
    ) where {FT}

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
update!(air::AirLayer{FT};
        f_CO₂::Union{Number,Nothing} = nothing,
        p_CO₂::Union{Number,Nothing} = nothing,
        p_H₂O::Union{Number,Nothing} = nothing,
        rh::Union{Number,Nothing} = nothing,
        t::Union{Number,Nothing} = nothing,
        vpd::Union{Number,Nothing} = nothing,
        wind::Union{Number,Nothing} = nothing
) where {FT} = (
    if !isnothing(t) air.auxil.t = t; end;
    if !isnothing(wind) air.auxil.wind = wind; end;
    if !isnothing(f_CO₂) air.auxil.f_CO₂ = f_CO₂; air.auxil.ps[2] = air.auxil.f_CO₂ * air.state.p_air * 1e-6; end;
    if !isnothing(p_CO₂) air.auxil.ps[2] = p_CO₂; air.auxil.f_CO₂ = air.auxil.ps[2] / air.state.p_air * 1e6; end;
    if !isnothing(p_H₂O) air.auxil.ps[3] = p_H₂O; end;
    if !isnothing(rh) air.auxil.ps[3] = saturation_vapor_pressure(air.auxil.t) * rh; end;
    if !isnothing(vpd) air.auxil.ps[3] = max(0, saturation_vapor_pressure(air.auxil.t) - vpd); end;

    return nothing
);


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
    (; AIRS, BRANCHES, CANOPY, LEAVES, ROOTS, SOIL_BULK, SOILS, TRUNK) = spac;

    # update chlorophyll and carotenoid contents (and spectra)
    if !isnothing(cab)
        for leaf in LEAVES
            leaf.bio.state.cab = cab;
        end;
    end;
    if !isnothing(car)
        for leaf in LEAVES
            leaf.bio.state.car = car;
        end;
    end;
    if !isnothing(cab) || !isnothing(car)
        plant_leaf_spectra!(config, spac);
    end;

    # update LAI and Vcmax (with scaling factor)
    if !isnothing(lai)
        CANOPY.structure.state.lai = lai;
        CANOPY.structure.state.δlai = lai .* ones(FT, DIM_LAYER) ./ DIM_LAYER;
        CANOPY.structure.auxil.x_bnds = (lai == 0 ? (collect(0:DIM_LAYER) ./ -DIM_LAYER) : ([0; [sum(CANOPY.structure.state.δlai[1:i]) for i in 1:DIM_LAYER]] ./ -lai));
        for i in 1:DIM_LAYER
            LEAVES[i].xylem.state.area = SOIL_BULK.state.area * CANOPY.structure.state.δlai[i];
        end;

        # make sure leaf area index setup and energy are correct
        for i in eachindex(LEAVES)
            LEAVES[i].xylem.state.area = SOIL_BULK.state.area * CANOPY.structure.state.δlai[i];
            initialize_struct!(LEAVES[i]);
        end;
    end;
    if !isnothing(vcmax)
        LEAVES[1].photosystem.state.v_cmax25 = vcmax;
    end;
    if !isnothing(vcmax) || !isnothing(lai)
        for i in 2:DIM_LAYER
            _scaling = isnothing(vcmax_expo) ? 1 : exp(-vcmax_expo * sum(CANOPY.structure.state.δlai[1:i-1]));
            LEAVES[i].photosystem.state.v_cmax25 = LEAVES[1].photosystem.state.v_cmax25 * _scaling;
            LEAVES[i].photosystem.state.j_max25 = LEAVES[1].photosystem.state.v_cmax25 * 1.67 * _scaling;
            LEAVES[i].photosystem.state.r_d25 = LEAVES[1].photosystem.state.v_cmax25 * 0.015 * _scaling;
            LEAVES[i].photosystem.auxil._t = 0;
        end;
    end;

    # update CI
    if !isnothing(ci)
        CANOPY.structure.auxil.ci = ci;
        CANOPY.structure.state.Ω_A = ci;
        CANOPY.structure.state.Ω_B = 0;
    end;

    # update Vcmax and Jmax TD
    if !isnothing(t_clm)
        for leaf in LEAVES
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
            _trunk_percent = TRUNK.HS.ΔH / (TRUNK.HS.ΔH + BRANCHES[end].HS.ΔH);
            (2 * kmax, 4 * kmax / _trunk_percent, 4 * kmax / (1 - _trunk_percent), 4 * kmax)
        else
            @assert length(kmax) == 4 "kmax must be a number or a tuple of length 4";
            kmax
        end;

        # partition kmax into the roots based on xylem area
        for _ilayer in ROOTS
            _ilayer.HS.K_X = _ilayer.HS.AREA / TRUNK.HS.AREA * _ks[1] * _ilayer.HS.L / _ilayer.HS.AREA;
        end;
        TRUNK.HS.K_X = _ks[2] * TRUNK.HS.L / TRUNK.HS.AREA;
        for _ilayer in BRANCHES
            _ilayer.HS.K_X = _ilayer.HS.AREA / TRUNK.HS.AREA * _ks[3] * _ilayer.HS.L / _ilayer.HS.AREA;
        end;
        for _ilayer in LEAVES
            _ilayer.HS.K_SLA = _ks[4] / (CANOPY.structure.state.lai * SOIL_BULK.state.area);
        end;
    end;

    # prescribe soil water content (within [Θ_RES,Θ_SAT])
    if !isnothing(swcs)
        for i in eachindex(swcs)
            SOILS[i].state.θ = max(SOILS[i].state.vc.Θ_RES + eps(FT), min(SOILS[i].state.vc.Θ_SAT - eps(FT), swcs[i]));
            initialize_struct!(SOILS[i], AIRS[1]);
        end;
    end;

    # prescribe soil temperature
    if !isnothing(t_soils)
        for i in eachindex(t_soils)
            SOILS[i].auxil.t = t_soils[i];
            initialize_struct!(SOILS[i], AIRS[1]);
        end;
    end;

    # prescribe leaf temperature
    if !isnothing(t_leaf)
        for leaf in LEAVES
            leaf.energy.auxil.t = t_leaf;
            leaf.energy.auxil.cp = heat_capacitance(leaf);
            leaf.energy.state.Σe = leaf.energy.auxil.cp * leaf.energy.auxil.t;
        end;

        # make sure leaf area index setup and energy are correct
        for i in eachindex(LEAVES)
            LEAVES[i].xylem.state.area = SOIL_BULK.state.area * CANOPY.structure.state.δlai[i];
            initialize_struct!(LEAVES[i]);
        end;
    end;

    return nothing
);
