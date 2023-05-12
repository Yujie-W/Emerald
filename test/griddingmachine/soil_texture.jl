using JSON

using Emerald.EmeraldData.GriddingMachineData: deploy_griddingmachine_artifacts!, griddingmachine_json!, reprocess_data!


# GLDAS2 soil texture map
for NX in ["1X", "4X"]
    @show NX;
    texture_json = "$(@__DIR__)/../../json/SOIL_TEXTURE_$(NX)_1Y_V1.json";

    if !isfile(texture_json)
        griddingmachine_json!(texture_json);
    end;

    json_dict = JSON.parse(open(texture_json));
end;
