
import Emerald.EmeraldData.GlobalDatasets as GD
import Emerald.EmeraldData.WeatherDrivers as WD
import Emerald.EmeraldFrontier as EF
import Emerald.EmeraldLand.Namespace as NS
import Emerald.EmeraldLand.SPAC as SPAC


begin "QL Model"

    gm_tag = "gm2";
    wd_tag = "wd1";
    gm_dict = GD.grid_dict(GD.LandDatasetLabels(gm_tag, 2019), 38.74, -92.20);
    gm_dict["MESSAGE_LEVEL"] = 2;
    wdf = WD.grid_weather_driver("wd1", gm_dict);
    for label in EF.DF_VARIABLES
        wdf[!,label] .= 0.0;
    end;
    for label in EF.DF_SIMULATIONS
        wdf[!,label] .= NaN;
    end;
    wdf_q = deepcopy(wdf);

    config = EF.spac_config(gm_dict);
    spac_q = GD.grid_spac(config, gm_dict);
    for r in [spac_q.plant.roots; spac_q.plant.trunk; spac_q.plant.branches; spac_q.plant.leaves]
        r.xylem.trait.vc.B = 5.703;
        r.xylem.trait.vc.C = 0.953;
    end;
    for l in spac_q.plant.leaves
        l.flux.trait.stomatal_model = NS.WangSM{gm_dict["FT"]}();
        l.photosystem.trait.FLM = NS.QLFluorescenceModelHanC3(gm_dict["FT"]);
    end;
    SPAC.prescribe_traits!(config, spac_q; kmax = 0.4);

    EF.simulation!(config, spac_q, wdf; selection = 3673:3674);

end;
