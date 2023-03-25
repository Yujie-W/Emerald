#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-25: add function to try the
#
#######################################################################################################################################################################################################
"""

    prescribe!(spac::MultiLayerSPAC{FT}, dfr::DataFrameRow) where {FT<:AbstractFloat}

Prescribe traits and environmental conditions, given
- `spac` `MultiLayerSPAC` type SPAC
- `dfr` `DataFrameRow` type weather driver

"""
function prescribe!(spac::MultiLayerSPAC{FT}, dfr::DataFrameRow) where {FT<:AbstractFloat}
    # read the data out of dataframe row to reduce memory allocation
    _df_atm::FT = dfr.P_ATM;
    _df_chl::FT = dfr.Chlorophyll;
    _df_cli::FT = dfr.CI;
    _df_co2::FT = dfr.CO2;
    _df_dif::FT = dfr.RAD_DIF;
    _df_dir::FT = dfr.RAD_DIR;
    _df_doy::FT = dfr.FDOY;
    _df_lai::FT = dfr.LAI;
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
    _df_vcm::FT = dfr.Vcmax;
    _df_vpd::FT = dfr.VPD;
    _df_wnd::FT = dfr.WIND;

    # adjust optimum t based on 10 day moving average skin temperature, prescribe soil water contents and leaf temperature (for version B1 only)
    push!(spac.MEMORY.tem, _df_tlf);
    if length(spac.MEMORY.tem) > 240 deleteat!(spac.MEMORY.tem,1) end;
    update!(spac; swcs = (_df_sw1, _df_sw2, _df_sw3, _df_sw4), t_clm = nanmean(spac.MEMORY.tem), t_leaf = max(_df_tar, _df_tlf), t_soils = (_df_ts1, _df_ts2, _df_ts3, _df_ts4));

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
#     2023-Mar-25: add function to try the
#
#######################################################################################################################################################################################################
function simulation! end

simulation!(wd_tag::String, gm_dict::Dict{String,Any}; appending::Bool = false, displaying::Bool = false) = (
    _spac = spac(gm_dict);
    _wd = weather_driver(wd_tag, gm_dict; appending = appending, displaying = displaying);

    for _i in 1:10
        soil_plant_air_continuum!(_spac, 360; p_on = false, t_on = false, θ_on = false);
    end;

    return nothing
);
