#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-25: add function to prescribe parameters from weather drivers
#     2023-Mar-27: prescribe T only if t_on is true, prescribe SWC only is θ_on is true
#
#######################################################################################################################################################################################################
"""

    prescribe!(spac::MultiLayerSPAC{FT}, dfr::DataFrameRow; t_on::Bool = true, θ_on::Bool = true) where {FT<:AbstractFloat}

Prescribe traits and environmental conditions, given
- `spac` `MultiLayerSPAC` type SPAC
- `dfr` `DataFrameRow` type weather driver
- `t_on` If true, plant energy budget is on, do not prescribe the temperatures
- `θ_on` If true, soil water budget is on, do not prescribe soil water contents

"""
function prescribe!(spac::MultiLayerSPAC{FT}, dfr::DataFrameRow; t_on::Bool = true, θ_on::Bool = true) where {FT<:AbstractFloat}
    # read the data out of dataframe row to reduce memory allocation
    _df_atm::FT = dfr.P_ATM;
    _df_chl::FT = dfr.CHLOROPHYLL;
    _df_cli::FT = dfr.CI;
    _df_co2::FT = dfr.CO2;
    _df_dif::FT = dfr.RAD_DIF;
    _df_dir::FT = dfr.RAD_DIR;
    _df_doy::FT = dfr.FDOY;
    _df_lai::FT = dfr.LAI;
    _df_pcp::FT = dfr.PRECIP;
    _df_sw1::FT = dfr.SWC_1;
    _df_sw2::FT = dfr.SWC_2;
    _df_sw3::FT = dfr.SWC_3;
    _df_sw4::FT = dfr.SWC_4;
    _df_tar::FT = dfr.T_AIR;
    _df_tlf::FT = dfr.T_LEAF;
    _df_ts1::FT = dfr.T_SOIL_1;
    _df_ts2::FT = dfr.T_SOIL_2;
    _df_ts3::FT = dfr.T_SOIL_3;
    _df_ts4::FT = dfr.T_SOIL_4;
    _df_vcm::FT = dfr.VCMAX25;
    _df_vpd::FT = dfr.VPD;
    _df_wnd::FT = dfr.WIND;

    # adjust optimum t based on 10 day moving average skin temperature
    _tleaf = t_on ? nanmean([_layer.t for _layer in spac.LEAVES]) : _df_tlf;
    push!(spac.MEMORY.tem, _tleaf);
    if length(spac.MEMORY.tem) > 240 deleteat!(spac.MEMORY.tem,1) end;
    update!(spac; t_clm = nanmean(spac.MEMORY.tem));

    # prescribe soil water contents and leaf temperature (for version B1 only)
    if !t_on
        update!(spac; t_leaf = max(_df_tar, _df_tlf), t_soils = (_df_ts1, _df_ts2, _df_ts3, _df_ts4));
    end;
    if !θ_on
        update!(spac; swcs = (_df_sw1, _df_sw2, _df_sw3, _df_sw4));
    end;

    # prescribe the precipitation related parameters
    spac.METEO.rain = _df_pcp * ρ_H₂O() / M_H₂O() / 3600;
    spac.METEO.t_precip = _df_tar;

    # if total LAI, Vcmax, or Chl changes, update them (add vertical Vcmax profile as well)
    _trigger_lai::Bool = !isnan(_df_lai) && (_df_lai != spac.MEMORY.lai);
    _trigger_vcm::Bool = !isnan(_df_vcm) && (_df_vcm != spac.MEMORY.vcm);
    _trigger_chl::Bool = !isnan(_df_chl) && (_df_chl != spac.MEMORY.chl);
    if _trigger_lai
        update!(spac; lai = _df_lai, vcmax_expo = 0.3);
        spac.MEMORY.lai = _df_lai;
    end;

    if _trigger_vcm
        update!(spac; vcmax = _df_vcm, vcmax_expo = 0.3);
        spac.MEMORY.vcm = _df_vcm;
    end;

    if _trigger_chl
        update!(spac; cab = _df_chl, car = _df_chl / 7);
        spac.MEMORY.chl = _df_chl;
    end;

    # update clumping index
    spac.CANOPY.ci = _df_cli;
    spac.CANOPY.Ω_A = _df_cli;

    # update environmental conditions
    for _alayer in spac.AIR
        _alayer.P_AIR = _df_atm;
        update!(_alayer; f_CO₂ = _df_co2, t = _df_tar, vpd = _df_vpd, wind = _df_wnd);
    end;

    # update shortwave radiation
    _in_dir = spac.RAD_SW_REF.e_direct' * spac.CANOPY.WLSET.ΔΛ / 1000;
    _in_dif = spac.RAD_SW_REF.e_diffuse' * spac.CANOPY.WLSET.ΔΛ / 1000;
    spac.RAD_SW.e_direct  .= spac.RAD_SW_REF.e_direct  .* _df_dir ./ _in_dir;
    spac.RAD_SW.e_diffuse .= spac.RAD_SW_REF.e_diffuse .* _df_dif ./ _in_dif;

    # update solar zenith angle based on the time
    spac.ANGLES.sza = solar_zenith_angle(spac.LATITUDE, FT(_df_doy));

    return nothing
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-25: move function from ClimaLand-0.2
#     2023-Mar-25: set reflectance based value to NaN at night
#     2023-Mar-27: add p_on, t_on, and θ_on options as in spac! function
#     2023-Mar-28: if option saving is false, return the simulated result dataframe
#     2023-Mar-28: add option selection to run part of the whole year simulations
#     2023-Mar-28: save swcs and temperatures based on t_on and θ_on
#
#######################################################################################################################################################################################################
"""

    simulation!(wd_tag::String,
                gm_dict::Dict{String,Any};
                appending::Bool = false,
                displaying::Bool = false,
                p_on::Bool = true,
                saving::Union{Nothing,String} = nothing,
                selection = :,
                t_on::Bool = true,
                θ_on::Bool = true)

Run simulation on site level, given
- `wd_tag` Weather drive tag such as `wd1`
- `gm_dict` GriddingMachine dict for site information
- `appending` If true, append new variables to weather driver when querying the file (set it to true when encountering any errors)
- `displaying` If true, displaying information regarding the steps
- `p_on` If true, plant hydraulic flow and pressure profiles will be updated
- `saving` If is not nothing, save the simulations as a Netcdf file in the working directory; if is nothing, return the simulated result dataframe
- `selection` Run selection of data, default is : (namely 1:end)
- `t_on` If true, plant energy budget is on, do not prescribe the temperatures
- `θ_on` If true, soil water budget is on, do not prescribe soil water contents

"""
function simulation! end

simulation!(wd_tag::String,
            gm_dict::Dict{String,Any};
            appending::Bool = false,
            displaying::Bool = false,
            p_on::Bool = true,
            saving::Union{Nothing,String} = nothing,
            selection = :,
            t_on::Bool = true,
            θ_on::Bool = true) = (
    _spac = spac(gm_dict);
    _wdf = weather_driver(wd_tag, gm_dict; appending = appending, displaying = displaying);
    _wdfr = eachrow(_wdf);

    # iterate through the time steps
    @showprogress for _dfr in _wdfr[selection]
        simulation!(_spac, _dfr; p_on = p_on, t_on = t_on, θ_on = θ_on);
    end;

    # save simulation results to hard drive
    if !isnothing(saving)
        save_nc!(saving, _wdf[selection,:]);

        return nothing
    end;

    return _wdf
);

simulation!(spac::MultiLayerSPAC{FT}, dfr::DataFrameRow; n_step::Int = 10, p_on::Bool = true, t_on::Bool = true, δt::Number = 3600, θ_on::Bool = true) where {FT<:AbstractFloat} = (
    # read the data out of dataframe row to reduce memory allocation
    _df_dif::FT = dfr.RAD_DIF;
    _df_dir::FT = dfr.RAD_DIR;

    # prescribe parameters
    prescribe!(spac, dfr; t_on = t_on, θ_on = θ_on);

    # run the model
    for _ in 1:n_step
        soil_plant_air_continuum!(spac, δt / n_step; p_on = p_on, t_on = t_on, θ_on = θ_on);
    end;

    # save the SIF and reflectance if there is sunlight
    if _df_dir + _df_dif >= 10
        dfr.BLUE   = MODIS_BLUE(spac.CANOPY);
        dfr.EVI    = MODIS_EVI(spac.CANOPY);
        dfr.NDVI   = MODIS_NDVI(spac.CANOPY);
        dfr.NIR    = MODIS_NIR(spac.CANOPY);
        dfr.NIRvI  = MODIS_NIRv(spac.CANOPY);
        dfr.NIRvR  = MODIS_NIRvR(spac.CANOPY);
        dfr.PAR    = spac.CANOPY.RADIATION.par_in;
        dfr.PPAR   = PPAR(spac);
        dfr.RED    = MODIS_RED(spac.CANOPY);
        dfr.SIF683 = TROPOMI_SIF683(spac.CANOPY);
        dfr.SIF740 = TROPOMI_SIF740(spac.CANOPY);
        dfr.SIF757 = OCO2_SIF759(spac.CANOPY);
        dfr.SIF771 = OCO2_SIF770(spac.CANOPY);
    else
        dfr.BLUE   = NaN;
        dfr.EVI    = NaN;
        dfr.NDVI   = NaN;
        dfr.NIR    = NaN;
        dfr.NIRvI  = NaN;
        dfr.NIRvR  = NaN;
        dfr.RED    = NaN;
    end;

    # save water contents and temperatures based on t_on and θ_on
    if θ_on
        dfr.MOD_SWC_1 = spac.SOIL.LAYERS[1].θ;
        dfr.MOD_SWC_2 = spac.SOIL.LAYERS[2].θ;
        dfr.MOD_SWC_3 = spac.SOIL.LAYERS[3].θ;
        dfr.MOD_SWC_4 = spac.SOIL.LAYERS[4].θ;
    end;
    if t_on
        _tleaf = [_leaf.t for _leaf in spac.LEAVES];
        dfr.MOD_T_L_MAX  = nanmax(_tleaf);
        dfr.MOD_T_L_MEAN = nanmean(_tleaf);
        dfr.MOD_T_L_MIN  = nanmin(_tleaf);
        dfr.MOD_T_S_1    = spac.SOIL.LAYERS[1].t;
        dfr.MOD_T_S_2    = spac.SOIL.LAYERS[2].t;
        dfr.MOD_T_S_3    = spac.SOIL.LAYERS[3].t;
        dfr.MOD_T_S_4    = spac.SOIL.LAYERS[4].t;
    end;

    # save the total flux into the DataFrame
    dfr.BETA  = BETA(spac);
    dfr.F_CO2 = CNPP(spac);
    dfr.F_GPP = GPP(spac);
    dfr.F_H2O = T_VEG(spac);

    return nothing
);
