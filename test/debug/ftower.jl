using DataFrames: DataFrame
using PkgUtility.DataIO: read_csv
using GriddingMachine.Indexer: CO₂_ppm

using Emerald.Land: GPP, LONGWAVE_OUT, SHORTWAVE_OUT
using Emerald.Land: prescribe!
using Emerald.Namespace: BulkSPAC, SPACConfig
using Emerald.SPAC: initialize_spac!, prescribe_traits!, push_t_history!, soil_plant_air_continuum!


# prepare the weather drivers from the csv files
function prepare_wd_nt(csvfile::String, utc_offset::Number)
    # things to prepare include the following:
    #     - driver_doy::FT = driver.FDOY[ind];
    #     - driver_wnd::FT = driver.WIND[ind];
    #     - driver_tar::FT = driver.TAIR[ind];
    #     - driver_vpd::FT = driver.VPD[ind];
    #     - driver_dif::FT = driver.RAD_SW_DIF[ind];
    #     - driver_dir::FT = driver.RAD_SW_DIR[ind];
    #     - driver_pcp::FT = driver.PPT[ind];
    #     - driver_chl::FT = driver.CHL[ind];
    #     - driver_lai::FT = driver.LAI[ind];
    #     - driver_cli::FT = driver.CI[ind];
    #     - driver_co2::FT = driver.CO2[ind];
    # currently not yet in the csv file
    #     - driver_atm::FT = driver.PATM[ind];
    #     - driver_lwr::FT = driver.RAD_LW[ind];
    #     - driver_b6f::FT = driver.B6F[ind];
    #     - driver_jmx::FT = driver.JMAX25[ind];
    #     - driver_vcm::FT = driver.VCMAX25[ind];
    indf = read_csv(csvfile);
    newdf = DataFrame();
    newdf[!,"FDOY"]       = @. indf.DOY + (indf.hour + 0.25 + indf.longitude / 15 - utc_offset) / 24;
    newdf[!,"WIND"]       = indf.Wind;
    newdf[!,"TAIR"]       = @. indf.Tair_C + 273.15;
    newdf[!,"VPD"]        = @. indf.VPD_kPa * 1000;
    newdf[!,"RAD_SW_DIR"] = @. indf.SW_IN * 0.8987696825369167;
    newdf[!,"RAD_SW_DIF"] = @. indf.SW_IN * 0.1012303174630833;
    newdf[!,"PPT"]        = @. indf.Precip / 1000;
    newdf[!,"LAI"]        = indf.LAI;
    newdf[!,"CHL"]        = indf.LCC;
    newdf[!,"CI"]         = indf.CI;
    newdf[!,"CO2"]        = ones(Float64, length(indf.DOY));
    # TODO
    newdf[!,"PATM"]       = ones(Float64, length(indf.DOY)) .* 101325;
    newdf[!,"RAD_LW"]     = ones(Float64, length(indf.DOY)) .* 280;
    newdf[!,"VCMAX25"]    = ones(Float64, length(indf.DOY)) .* 50;
    newdf[!,"JMAX25"]     = ones(Float64, length(indf.DOY)) .* 100;
    newdf[!,"B6F"]        = ones(Float64, length(indf.DOY)) .* 0.5;

    # update CO₂ based on the year
    for yy in unique(indf.year)
        newdf.CO2[indf.year .== yy] .= CO₂_ppm(yy);
    end;

    # make sure there is no NaN in the data
    for k in names(newdf)
        if any(isnan.(newdf[!,k]))
            @warn "NaN values found in column $(k) of $(csvfile)...";
        end;
    end;

    return NamedTuple{Tuple(Symbol.(names(newdf)))}(Tuple([newdf[!,k] for k in names(newdf)]))
end;


# nt = prepare_wd_nt("/mnt/c/Users/wyujie/Desktop/EC_sunsh_input/AR-SLu_input.csv", -3);
# wd = GriddingMachine.Indexer.grid_weather("wd1", 2009, -66.4598, -33.4648);


function run_model_remi!(csvfile::String, lat::FT, lon::FT, δt::Number, utc_offset::Number) where {FT<:AbstractFloat}
    driver = prepare_wd_nt(csvfile, utc_offset);

    #
    # config
    #
    CONFIG = SPACConfig(FT);
    CONFIG.FEATURES.ENABLE_REF = false;  # disable reflectance spectrum
    CONFIG.FEATURES.ENABLE_SIF = false;  # disable solar-induced fluorescence
    CONFIG.CONFIG_INFO.MESSAGE_LEVEL = 2;
    CONFIG.FEATURES.UNLIMITED_NSC_POOL = true;

    #
    # spac
    #
    spac = BulkSPAC(CONFIG; latitude = lat, longitude = lon);
    spac.soil_bulk.trait.color = 1;
    spac.plant.pool.c_pool = Inf;
    prescribe_traits!(CONFIG, spac; sai = 0);

    prescribe!(CONFIG, spac, driver, 1; initialize_state = true);
    for ind in eachindex(driver.FDOY)[1:8760]
        @show ind;
        prescribe!(CONFIG, spac, driver, ind);
        soil_plant_air_continuum!(CONFIG, spac, δt);
        push_t_history!(CONFIG, spac);
    end;

    return nothing
end;


run_model_remi!("/mnt/c/Users/wyujie/Desktop/EC_sunsh_input/AR-SLu_input.csv", -66.4598, -33.4648, 1800, -3);
