using JSON

using Emerald.EmeraldData.GriddingMachineData: deploy_griddingmachine_artifacts!, griddingmachine_json!, reprocess_data!


# TROPOMI SIF at 740 nm
for NX in ["1X", "2X", "4X", "5X", "12X"]
    for MT in ["1M", "8D"]
        for DC in ["", "DC"]
            tropomi_json = "$(@__DIR__)/../../json/TROPOMI_740$(DC)_$(NX)_$(MT)_V1.json";

            if !isfile(tropomi_json)
                griddingmachine_json!(tropomi_json);
            end;

            json_dict = JSON.parse(open(tropomi_json));
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
end;
