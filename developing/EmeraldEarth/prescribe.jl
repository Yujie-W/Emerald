#=
# this function is a tailored function meant to create a SPAC for a given date
site_spac(config::SPACConfig{FT}, gmd::Dict{String,Any}, iday::Int) where {FT} = (
    spac = site_spac(config, gmd);
    prescribe_gm_wd_data!(config, spac, gmd, iday);

    return spac
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the spac
#     2023-Mar-13: add step to synchronize state variables into spac
#     2023-Mar-29: prescribe longwave radiation as well
#     2023-Jun-15: make sure prescribed swc does not exceed the limits
#     2023-Jun-15: make sure prescribed soil parameters are not NaN and rad is >= 0
#     2024-Feb-28: rename function to prescribe_gm_wd_data! and move it to EmeraldData
#     2024-Feb-28: make it possible to initialize spac with ss_dict
#     2024-Apr-17: update solar azimuth angle as well
#     2025-Mar-13: add method to update SPAC for a given date (without weather driver data)
#
#######################################################################################################################################################################################################
"""

    prescribe_gm_wd_data!(config::SPACConfig{FT}, spac::BulkSPAC{FT}, gmd::Dict{String,Any}, wd_dict::Dict{String,Any}, ss_dict::Union{Dict{String,Any},Nothing} = nothing) where {FT}

Prescribe the SPAC with GriddingMachine and weather driver data, given
- `config` Configurations for SPAC
- `spac` SPAC to be prescribed
- `gmd` Dictionary of GriddingMachine data
- `wd_dict` Dictionary of weather driver data
- `ss_dict` Dictionary of initial states

"""
function prescribe_gm_wd_data! end;

prescribe_gm_wd_data!(config::SPACConfig{FT}, spac::BulkSPAC{FT}, gmd::Dict{String,Any}, wd_dict::Dict{String,Any}, ss_dict::Union{Dict{String,Any},Nothing} = nothing) where {FT} = (
    # update environmental conditions
    for air in spac.airs
        air.state.p_air = wd_dict["P_ATM"];
        prescribe_air!(air; f_CO₂ = gmd["CO2"], t = wd_dict["T_AIR"], vpd = wd_dict["VPD"], wind = wd_dict["WIND"]);
    end;

    # update shortwave and longwave radiation
    ref_dir = view(config.CONSTANTS.SPECTRA.SOLAR_RAD,:,1)'  * config.CONSTANTS.SPECTRA.ΔΛ / 1000;
    ref_dif = view(config.CONSTANTS.SPECTRA.SOLAR_RAD,:,2)' * config.CONSTANTS.SPECTRA.ΔΛ / 1000;
    spac.meteo.rad_sw.e_dir .= view(config.CONSTANTS.SPECTRA.SOLAR_RAD,:,1) .* max(0,wd_dict["RAD_DIR"]) ./ ref_dir;
    spac.meteo.rad_sw.e_dif .= view(config.CONSTANTS.SPECTRA.SOLAR_RAD,:,2) .* max(0,wd_dict["RAD_DIF"]) ./ ref_dif;
    spac.meteo.rad_lw = wd_dict["RAD_LW"];
    saa = solar_azimuth_angle(spac.info.lat, FT(wd_dict["FDOY"]));
    sza = solar_zenith_angle(spac.info.lat, FT(wd_dict["FDOY"]));
    spac.canopy.sun_geometry.state.saa = saa;
    spac.canopy.sun_geometry.state.sza = (wd_dict["RAD_DIR"] + wd_dict["RAD_DIF"] > 10) ? min(sza, 88) : sza;

    # update t_clm to make Vcmax25 and Jmax25 TD temperature dependent
    prescribe_traits!(config, spac; t_clm = mean(spac.plant.memory.t_history));

    # synchronize LAI, CHL, and CI
    iday = Int(floor(wd_dict["INDEX"] / 24)) + 1;
    chl = gmd["CHLOROPHYLL"][iday];
    ci = gmd["CLUMPING"][iday];
    lai = gmd["LAI"][iday];
    vcm = gmd["VCMAX25"][iday];
    prescribe_traits!(config, spac; cab = chl, car = chl / 7, ci = ci, lai = lai, vcmax = vcm, vertical_expo = 0.3);

    # if ss_dict is not nothing, update soil water content and leaf temperature
    if !isnothing(ss_dict)
        # update soil water content
        swckeys = ["SWC_1", "SWC_2", "SWC_3", "SWC_4"];
        tslkeys = ["T_S_1", "T_S_2", "T_S_3", "T_S_4"];
        prescribe_soil!(spac; swcs = Tuple(min(spac.soils[i].trait.vc.Θ_SAT - 0.001, ss_dict[swckeys[i]]) for i in 1:4), t_soils = Tuple(ss_dict[tslkeys[i]] for i in 1:4));

        # prescribe leaf temperature from skin temperature
        @. spac.plant.memory.t_history = ss_dict["T_SKN"];
        prescribe_traits!(config, spac; t_leaf = ss_dict["T_SKN"], t_clm = mean(spac.plant.memory.t_history));
    end;

    return nothing
);

# this function is a tailored function meant to create a SPAC for a given date
prescribe_gm_wd_data!(config::SPACConfig{FT}, spac::BulkSPAC{FT}, gmd::Dict{String,Any}, iday::Int) where {FT} = (
    # synchronize LAI, CHL, and CI
    chl = gmd["CHLOROPHYLL"][iday];
    ci = gmd["CLUMPING"][iday];
    lai = gmd["LAI"][iday];
    vcm = gmd["VCMAX25"][iday];
    prescribe_traits!(config, spac; cab = chl, car = chl / 7, ci = ci, lai = lai, vcmax = vcm, vertical_expo = 0.3);

    return nothing
);
=#
