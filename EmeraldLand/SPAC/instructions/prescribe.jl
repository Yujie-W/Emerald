# This file contains functions to prescribe environmental variables and traits

#######################################################################################################################################################################################################
#
# Changes to this method
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
