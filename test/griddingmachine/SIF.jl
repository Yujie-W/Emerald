#=
using JSON

using Emerald.EmeraldData.GriddingMachineData: deploy_griddingmachine_artifacts!, griddingmachine_json!, reprocess_data!


# TROPOMI SIF at 740 nm
for NX in [1, 2, 4, 5, 12]
    for MT in ["1M", "8D"]
        # for 740 (not DC)
        tropomi_json = "$(@__DIR__)/../../json/TROPOMI_740_V1.json";
        json_dict = JSON.parse(open(tropomi_json));
        json_dict["GRIDDINGMACHINE"]["LAT_LON_RESO"] = NX;
        json_dict["GRIDDINGMACHINE"]["TEMPORAL_RESO"] = MT;
        json_dict["INPUT_MAP_SETS"]["FILE_NAME_PATTERN"] = "SIF_TROPOMI_740_$(NX)X_$(MT)_XXXXXXXX_V1.nc";

        name_function = eval(Meta.parse(json_dict["INPUT_MAP_SETS"]["FILE_NAME_FUNCTION"]));
        data_scaling_functions = [_dict["SCALING_FUNCTION"] == "" ? nothing :  eval(Meta.parse(_dict["SCALING_FUNCTION"])) for _dict in json_dict["INPUT_VAR_SETS"]];
        std_scaling_functions = if "INPUT_STD_SETS" in keys(json_dict)
            [_dict["SCALING_FUNCTION"] == "" ? nothing : eval(Meta.parse(_dict["SCALING_FUNCTION"])) for _dict in json_dict["INPUT_STD_SETS"]]
        else
            [nothing for _dict in json_dict["INPUT_VAR_SETS"]]
        end;
        reprocess_data!(json_dict; file_name_function = name_function, data_scaling_functions = data_scaling_functions, std_scaling_functions = std_scaling_functions);
        deploy_griddingmachine_artifacts!(json_dict);

        # for 740 DC
        tropomi_json = "$(@__DIR__)/../../json/TROPOMI_740DC_V1.json";
        json_dict = JSON.parse(open(tropomi_json));
        json_dict["GRIDDINGMACHINE"]["LAT_LON_RESO"] = NX;
        json_dict["GRIDDINGMACHINE"]["TEMPORAL_RESO"] = MT;
        json_dict["INPUT_MAP_SETS"]["FILE_NAME_PATTERN"] = "SIF_TROPOMI_740_$(NX)X_$(MT)_XXXXXXXX_V1.nc";

        name_function = eval(Meta.parse(json_dict["INPUT_MAP_SETS"]["FILE_NAME_FUNCTION"]));
        data_scaling_functions = [_dict["SCALING_FUNCTION"] == "" ? nothing :  eval(Meta.parse(_dict["SCALING_FUNCTION"])) for _dict in json_dict["INPUT_VAR_SETS"]];
        std_scaling_functions = if "INPUT_STD_SETS" in keys(json_dict)
            [_dict["SCALING_FUNCTION"] == "" ? nothing : eval(Meta.parse(_dict["SCALING_FUNCTION"])) for _dict in json_dict["INPUT_STD_SETS"]]
        else
            [nothing for _dict in json_dict["INPUT_VAR_SETS"]]
        end;
        reprocess_data!(json_dict; file_name_function = name_function, data_scaling_functions = data_scaling_functions, std_scaling_functions = std_scaling_functions);
        deploy_griddingmachine_artifacts!(json_dict);
    end;
end;
=#
