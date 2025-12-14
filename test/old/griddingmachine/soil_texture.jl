#=
using JSON

using Emerald.EmeraldData.GriddingMachineData: deploy_griddingmachine_artifacts!, griddingmachine_json!, reprocess_data!


# GLDAS2 soil texture map
for NX in [1, 4]
    texture_json = "$(@__DIR__)/../../json/SOIL_TEXTURE_V1.json";
    json_dict = JSON.parse(open(texture_json));
    json_dict["GRIDDINGMACHINE"]["LAT_LON_RESO"] = NX;
    json_dict["INPUT_MAP_SETS"]["FILE_NAME_PATTERN"] = "SOIL_TEXTURE_$(NX)X_1Y_V1.nc4";

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
=#
