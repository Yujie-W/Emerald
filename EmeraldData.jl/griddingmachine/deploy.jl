"""

    deploy_griddingmachine_artifacts!(dict::Dict)

Deploy GriddingMachine artifacts, given
- `dict` A dictionary based on JSON file

"""
function deploy_griddingmachine_artifacts!(dict::Dict)
    _dict_grid = dict["GRIDDINGMACHINE"];

    # determine if there is any information for years
    _years = _dict_grid["YEARS"];
    _i_years = (_years == "" ? [1] : eachindex(_years));
    for _i_year in _i_years
        _tag = (_years == "" ? griddingmachine_tag(dict, 0) : griddingmachine_tag(dict, _years[_i_year]));
        _artifact_files = ["GRIDDINGMACHINE", "$(_tag).nc"];
        deploy_artifact!(ARTIFACT_TOML, _tag, DATASET_FOLDER, _artifact_files, ARTIFACT_FOLDER, FTP_URLS);
    end;

    return nothing
end
