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
function update! end


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
    ) where {FT<:AbstractFloat}

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
) where {FT<:AbstractFloat} = (
    if !isnothing(t) air.t = t; end;
    if !isnothing(wind) air.wind = wind; end;
    if !isnothing(f_CO₂) air.f_CO₂ = f_CO₂; air.p_CO₂ = air.f_CO₂ * air.P_AIR * 1e-6; end;
    if !isnothing(p_CO₂) air.p_CO₂ = p_CO₂; air.f_CO₂ = air.p_CO₂ / air.P_AIR * 1e6; end;
    if !isnothing(p_H₂O) air.p_H₂O = p_H₂O; end;
    if !isnothing(rh) air.p_H₂O = saturation_vapor_pressure(air.t) * rh; end;
    if !isnothing(vpd) air.p_H₂O = max(0, saturation_vapor_pressure(air.t) - vpd); end;

    return nothing
);


"""

    update!(spac::MultiLayerSPAC{FT},
            config::SPACConfiguration{FT};
            cab::Union{Number,Nothing} = nothing,
            car::Union{Number,Nothing} = nothing,
            ci::Union{Number,Nothing} = nothing,
            lai::Union{Number,Nothing} = nothing,
            swcs::Union{Tuple,Nothing} = nothing,
            t_clm::Union{Number,Nothing} = nothing,
            t_leaf::Union{Number,Nothing} = nothing,
            t_soils::Union{Tuple,Nothing} = nothing,
            vcmax::Union{Number,Nothing} = nothing,
            vcmax_expo::Union{Number,Nothing} = nothing
    ) where {FT<:AbstractFloat}

Update the physiological parameters of the SPAC, given
- `spac` Soil plant air continuum
- `config` Configuration for `MultiLayerSPAC`
- `cab` Chlorophyll content. Optional, default is nothing
- `car` Carotenoid content. Optional, default is nothing
- `ci` Clumping index. Optional, default is nothing
- `lai` Leaf area index. Optional, default is nothing
- `swcs` Soil water content at different layers. Optional, default is nothing
- `t_clm` Moving average temperature to update Vcmax and Jmax temperature dependencies. Optional, default is nothing
- `t_leaf` Leaf temperature. Optional, default is nothing
- `t_soils` Soil temperature at different layers. Optional, default is nothing
- `vcmax` Vcmax25 at the top of canopy. Optional, default is nothing
-`vcmax_expo` Exponential tuning factor to adjust Vcmax25. Optional, default is nothing

"""
update!(spac::MultiLayerSPAC{FT},
        config::SPACConfiguration{FT};
        cab::Union{Number,Nothing} = nothing,
        car::Union{Number,Nothing} = nothing,
        ci::Union{Number,Nothing} = nothing,
        lai::Union{Number,Nothing} = nothing,
        swcs::Union{Tuple,Nothing} = nothing,
        t_clm::Union{Number,Nothing} = nothing,
        t_leaf::Union{Number,Nothing} = nothing,
        t_soils::Union{Tuple,Nothing} = nothing,
        vcmax::Union{Number,Nothing} = nothing,
        vcmax_expo::Union{Number,Nothing} = nothing
) where {FT<:AbstractFloat} = (
    (; CANOPY, DIM_LAYER, LEAVES, SOIL) = spac;

    # update chlorophyll and carotenoid contents (and )
    if !isnothing(cab)
        for _leaf in LEAVES
            _leaf.BIO.cab = cab;
            _leaf.BIO._v_storage = 0;
        end;
    end;
    if !isnothing(car)
        for _leaf in LEAVES
            _leaf.BIO.car = car;
            _leaf.BIO._v_storage = 0;
        end;
    end;
    if !isnothing(cab) || !isnothing(car)
        leaf_spectra!(spac, config);
    end;

    # update LAI
    if !isnothing(lai)
        CANOPY.lai = lai;
        for _i in 1:DIM_LAYER
            LEAVES[_i].HS.AREA = SOIL.AREA * CANOPY.lai / DIM_LAYER;
        end;
    end;

    # update CI
    if !isnothing(ci)
        CANOPY.ci = ci;
        CANOPY.Ω_A = ci;
        CANOPY.Ω_B = 0;
    end;

    # prescribe soil water content
    if !isnothing(swcs)
        for _i in eachindex(swcs)
            _slayer = SOIL.LAYERS[_i];
            _slayer.θ = swcs[_i];
            _slayer.e = (_slayer.CP * _slayer.ρ + _slayer.θ * CP_L() * ρ_H₂O()) * _slayer.t;
        end;
    end;

    # prescribe soil temperature
    if !isnothing(t_soils)
        for _i in eachindex(t_soils)
            _slayer = SOIL.LAYERS[_i];
            _slayer.t = t_soils[_i];
            _slayer.e = (_slayer.CP * _slayer.ρ + _slayer.θ * CP_L() * ρ_H₂O()) * _slayer.t;
        end;
    end;

    # update Vcmax and Jmax TD
    if !isnothing(t_clm)
        for _leaf in LEAVES
            _leaf.PSM.TD_VCMAX.ΔSV = 668.39 - 1.07 * (t_clm - T₀());
            _leaf.PSM.TD_JMAX.ΔSV = 659.70 - 0.75 * (t_clm - T₀());
        end;
    end;

    # prescribe leaf temperature
    if !isnothing(t_leaf)
        for _leaf in LEAVES
            _leaf.t = t_leaf;
            _leaf.e = (_leaf.CP * _leaf.BIO.lma * 10 + _leaf.HS.v_storage * CP_L_MOL(FT)) * _leaf.t;
        end;
    end;

    # update Vcmax at the top layer
    if !isnothing(vcmax)
        LEAVES[1].PSM.v_cmax25 = vcmax;
    end;

    # update Vcmax profile if lai or vcmax is given
    if !isnothing(vcmax) || !isnothing(lai)
        for _i in 2:DIM_LAYER
            _scaling = isnothing(vcmax_expo) ? 1 : exp(-vcmax_expo * CANOPY.lai * ((_i - 1) / DIM_LAYER));
            LEAVES[_i].PSM.v_cmax25 = LEAVES[1].PSM.v_cmax25 * _scaling;
            LEAVES[_i].PSM.j_max25 = LEAVES[1].PSM.v_cmax25 * 1.67 * _scaling;
            LEAVES[_i].PSM.r_d25 = LEAVES[1].PSM.v_cmax25 * 0.015 * _scaling;
            LEAVES[_i].PSM._t = 0;
        end;
    end;

    return nothing
);