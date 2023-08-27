using Emerald.EmeraldFrontier: DF_SIMULATIONS, DF_VARIABLES, GriddingMachineLabels, gm_dict, simulation!, weather_driver
using Emerald.EmeraldLand.Constant: GRAVITY, M_H₂O, T₀
using Emerald.EmeraldLand.Namespace: MultiLayerSPAC, SPACConfiguration, VanGenuchten, WeibullVC
using Emerald.EmeraldLand.SPAC: initialize!, soil_plant_air_continuum!, update!
using Emerald.EmeraldVisualization: canvas, decorate!, save_canvas!

using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: read_LUT

using NetcdfIO: read_nc


CONFIG = SPACConfiguration{Float64}(DEBUG = true);
MERRAFOLDER = "/home/wyujie/ProjectData/SIF_MIP/MERRA";


# function to load data from GriddingMachine using the functions provided by EmeraldFrontier
function prepare_gm_dict(site::String = "NR1", year::Int = 2019)
    @assert site in ["NR1", "OBS", "xDG"];
    @assert 2000 <= year <= 2021;

    # read in lat and lon
    if site == "NR1"
        _lat = 40.03;
        _lon = -105.55;
        _ele = 3050;
        _lai = 4.0;
        _vcmax25 = 30.0;
        _soil_n = [1.35301 for _ in 1:4];
        _soil_a = [129.66 for _ in 1:4];
        _soil_r = [0.196 for _ in 1:4];
        _soil_s = [0.477 for _ in 1:4];
    end;

    # read data from GriddingMachine
    _dict = gm_dict(GriddingMachineLabels(year = year), _lat, _lon)
    _dict["ELEVATION"] = _ele;
    _dict["LAI"] = _lai;
    _dict["VCMAX25"] = _vcmax25;
    _dict["SOIL_N"] = _soil_n;
    _dict["SOIL_α"] = _soil_a;
    _dict["SOIL_ΘR"] = _soil_r;
    _dict["SOIL_ΘS"] = _soil_s;

    # push extra variables into dict
    push!(_dict, "SITE" => site);

    return _dict
end;


# function to prepare weather drivers
function prepare_weather_driver(gmdict::Dict, wd::String = "MERRA")
    # use ERA5 data from Emerald because the data provided is problematic
    if wd == "ERA5"
        _df = weather_driver("wd1", gmdict);
        for _var in [DF_SIMULATIONS; DF_VARIABLES]
            _df[!,_var] .= NaN;
        end;

        return _df
    end;

    # read in data provided by SIF MIP
    _site = gmdict["SITE"];
    _year = gmdict["YEAR"];
    if wd == "MERRA"
        _df = read_nc("$(MERRAFOLDER)/ABOVE_MET_MERRA_$(_site)_$(_year)_ver.2022.07.14.nc");
    end;

    # reformate the data names to match the ones used in EmeraldFrontier
    _df.T_AIR = _df.ta .+ T₀();
    _df.RAD_LW = _df.lw;
    _df.RAD = _df.sw;
    _df.RAD_DIR = _df.pardir ./ max.(1e-3,_df.par) .* _df.sw;
    _df.RAD_DIF = _df.RAD .- _df.RAD_DIR;
    _df.P_ATM = _df.ps .* 100;
    _df.WIND = sqrt.(_df.uwnd .^ 2 .+ _df.vwnd .^ 2);
    _df.VPD = _df.vpd .* 100;

    return _df
end;


# functions to create SPAC for four sites, prepare the spac only when it is the first year of the simulation
function prepare_spac(gmdict::Dict)
    FT = gmdict["FT"];
    _lat = gmdict["LATITUDE"];
    _lon = gmdict["LONGITUDE"];
    _ele = gmdict["ELEVATION"];

    # create spac for NR1
    if gmdict["SITE"] == "NR1"
        # cerate SPAC, initialize here to make sure LAI and Kmax match
        _spac = MultiLayerSPAC(CONFIG; air_bounds = collect(0:0.5:14), elevation = _ele, latitude = _lat, longitude = _lon, zs = [-1,6,13]);
        _spac.SOIL.AREA = 32.088;
        initialize!(_spac, CONFIG);
        update!(_spac, CONFIG; ci = 0.48, lai = 4, vcmax = 30);

        # change vulnerability curve
        _vc = WeibullVC{FT}(4.09, 5.82);
        for _ilayer in [_spac.ROOTS; _spac.TRUNK; _spac.BRANCHES; _spac.LEAVES]
            _ilayer.HS.VC = deepcopy(_vc);
        end;
    end;

    # set plant kmax to 0.05 as in the GMD paper
    update!(_spac, CONFIG; kmax = 0.1);

    # change soil kmax to a higher value
    for _slayer in _spac.SOIL.LAYERS
        _slayer.VC.K_MAX = 1e-4 / GRAVITY(FT) * 1e6 / M_H₂O(FT);
    end;

    # set the limits of leaf gs within (0.0001, 0.1)
    for _ilayer in _spac.LEAVES
        _ilayer.G_LIMITS .= [1e-4, 0.1];
    end;

    # initialize the spac
    initialize!(_spac, CONFIG);

    return _spac
end;


# function to run CliMA Land for two years
function run_experiment_1!(; force::Bool = false)
    _sites = ["NR1"];
    for _site in _sites
        # create the spac to use in all the simulations
        _gmdict = prepare_gm_dict(_site, 2017);
        _spac = prepare_spac(_gmdict);

        # run the model for multiple years
        for _year in [2018]
            _modfile = "$(@__DIR__)/output/ERA5_$(_site)_$(_year).nc";
            _figfile = "$(@__DIR__)/output/ERA5_$(_site)_$(_year).png";
            if !isfile(_modfile) || force
                _gmdict = prepare_gm_dict(_site, _year);
                _wdf = prepare_weather_driver(_gmdict, "ERA5");
                simulation!(CONFIG, _spac, _wdf; initialial_state = false, saving = _modfile);
            end;
            if (isfile(_modfile) && !isfile(_figfile)) || force
                _fig = canvas(_site; figsize = (12,15), nrow = 8, ncol = 4);
                _df = read_nc(_modfile);
                _vars = [DF_SIMULATIONS; DF_VARIABLES];
                for _i in eachindex(_vars)
                    _fig.axes[_i].plot(_df.FDOY, _df[:,_vars[_i]]);
                end;
                decorate!(_fig.axes[1:length(_vars)]; add_title = false, yaxis_labels = _vars);
                save_canvas!(_fig, "ERA5_$(_site)_$(_year)"; folder = "$(@__DIR__)/output", formats = ["png"]);
            end;
        end;
    end;

    return nothing
end


# this part is to test why the simulation is not working
#=

site = "NR1";
# create the spac to use in all the simulations
gmdict = prepare_gm_dict(site, 2017);
spac = prepare_spac(gmdict);

# run the model for multiple years
year = 2018;
modfile = "$(@__DIR__)/output/ERA5_$(site)_$(year).nc";
gmdict = prepare_gm_dict(site, year);
wdf = prepare_weather_driver(gmdict, "ERA5");
wdfr = eachrow(wdf);
simulation!(CONFIG, spac, wdf; initialial_state = false, saving = modfile);


# make a deep copy of the spac to debug
# spac_bak = deepcopy(spac);


# read back in the spac from the backup
spac = deepcopy(spac_bak);
for dfr in wdfr[3860:end]
    simulation!(spac, CONFIG, dfr; p_on = true, t_on = true, θ_on = true);
end;

=#
