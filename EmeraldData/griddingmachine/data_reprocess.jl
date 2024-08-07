#=
#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-May-12: move function from GriddingMachineDatasets
#     2023-May-17: add coverage support
#     2023-May-17: add step to visualize reprocessed data before saving
#
#######################################################################################################################################################################################################
"""

    reprocess_data!(
                dict::Dict;
                file_name_function::Union{Function,Nothing} = nothing,
                data_scaling_functions::Vector = [nothing for i in eachindex(dict["VARIABLE_SETTINGS"])],
                std_scaling_functions::Vector = [nothing for i in eachindex(dict["VARIABLE_SETTINGS"])])

Reprocess the data to use in GriddingMachine artifacts, given
- `dict` JSON dict
- `file_name_function` Function to find file
- `data_scaling_functions` Functions to scale data
- `std_scaling_functions` Functions to scale std

"""
function reprocess_data!(
            dict::Dict;
            file_name_function::Union{Function,Nothing} = nothing,
            data_scaling_functions::Vector = [nothing for i in eachindex(dict["VARIABLE_SETTINGS"])],
            std_scaling_functions::Vector = [nothing for i in eachindex(dict["VARIABLE_SETTINGS"])])
    _dict_file = dict["INPUT_MAP_SETS"];
    _dict_grid = dict["GRIDDINGMACHINE"];
    _dict_vars = dict["INPUT_VAR_SETS"];
    _dict_outv = dict["OUTPUT_VAR_ATTR"];
    _dict_refs = dict["OUTPUT_REF_ATTR"];
    _dict_stds = "INPUT_STD_SETS" in keys(dict) ? dict["INPUT_STD_SETS"] : nothing;

    # determine if there is any information for years
    _years = _dict_grid["YEARS"];
    _files = [];
    if isnothing(_years)
        push!(_files, _dict_file["FOLDER"] * "/" * _dict_file["FILE_NAME_PATTERN"]);
    else
        for _year in _years
            push!(_files, _dict_file["FOLDER"] * "/" * replace(_dict_file["FILE_NAME_PATTERN"], "XXXXXXXX" => file_name_function(_year)));
        end;
    end;

    # iterate through the files
    _i_years = (isnothing(_years) ? [1] : eachindex(_years));
    for _i_year in _i_years
        # determine whether to skip based on the tag
        _tag = (isnothing(_years) ? griddingmachine_tag(dict) : griddingmachine_tag(dict, _years[_i_year]));
        _reprocessed_file = "/home/wyujie/GriddingMachine/reprocessed/$(_tag).nc";

        # reprocess the data only if file does not exist
        if !isfile(_reprocessed_file)
            @info "File $(_reprocessed_file) does not exist, reprocessing...";
            # read the data
            _file = _files[_i_year]
            if length(_dict_vars) == 1
                _reprocessed_data = read_data(
                            _file,
                            _dict_vars[1],
                            [_dict_file["FLIP_LAT"],_dict_file["FLIP_LON"]],
                            _dict_grid["LAT_LON_RESO"];
                            coverage = _dict_file["COVERAGE"],
                            scaling_function = data_scaling_functions[1]);
                _reprocessed_std = if !isnothing(_dict_stds)
                    read_data(
                            _file,
                            _dict_stds[1],
                            [_dict_file["FLIP_LAT"],_dict_file["FLIP_LON"]],
                            _dict_grid["LAT_LON_RESO"];
                            coverage = _dict_file["COVERAGE"],
                            scaling_function = std_scaling_functions[1])
                else
                    similar(_reprocessed_data) .* NaN
                end;
            else
                _reprocessed_data = ones(Float64, 360 * _dict_grid["LAT_LON_RESO"], 180 * _dict_grid["LAT_LON_RESO"], length(_dict_vars));
                _reprocessed_std = ones(Float64, 360 * _dict_grid["LAT_LON_RESO"], 180 * _dict_grid["LAT_LON_RESO"], length(_dict_vars)) .* NaN;
                for _i_var in eachindex(_dict_vars)
                    _reprocessed_data[:,:,_i_var] = read_data(
                            _file,
                            _dict_vars[_i_var],
                            [_dict_file["FLIP_LAT"],_dict_file["FLIP_LON"]],
                            _dict_grid["LAT_LON_RESO"];
                            coverage = _dict_file["COVERAGE"],
                            scaling_function = data_scaling_functions[_i_var]);
                    if !isnothing(_dict_stds)
                        _reprocessed_std[:,:,_i_var] = read_data(
                            _file,
                            _dict_stds[_i_var],
                            [_dict_file["FLIP_LAT"],_dict_file["FLIP_LON"]],
                            _dict_grid["LAT_LON_RESO"];
                            coverage = _dict_file["COVERAGE"],
                            scaling_function = std_scaling_functions[_i_var]);
                    end;
                end;
            end;

            # plot the reprocessed file to make sure it is processed correctly
            _reso = 1 / _dict_grid["LAT_LON_RESO"];
            _lats = collect((-90+_reso/2):_reso:90);
            _lons = collect((-180+_reso/2):_reso:180);
            _ffig = "$(@__DIR__)/temp.gif";
            # animate_data!(_lons, _lats, _reprocessed_data; filename = _ffig);
            _msg = "The figure is saved as $(_ffig). Is the generated data okay? Type Y/y or N/n to continue > ";
            _save_data = input_yes_or_no(_msg; bool_conversion = true);

            # save the file
            if _save_data
                _var_attr::Dict{String,String} = merge(_dict_outv,_dict_refs);
                _dim_names = length(size(_reprocessed_std)) == 3 ? ["lon", "lat", "ind"] : ["lon", "lat"];
                save_nc!(_reprocessed_file, "data", _reprocessed_data, _var_attr);
                append_nc!(_reprocessed_file, "std", _reprocessed_std, _var_attr, _dim_names);
            end;

            # delete the temporary figure
            rm(_ffig; force = true);
        else
            @info "File $(_reprocessed_file) exists already, skipping...";
        end;
    end;

    # add change logs based on the JSON file
    #=
    - Add uncertainty (filled with NaN)
    - Make the map to global scale (fill with NaN)
    - Reformatted from GeoTIFF or binary to NetCDF
    - Latitude and Longitude re-oriented to from South to North and from West to East
    - Data scaling removed (from log(x), exp(x), or kx+b to x)
    - Data regridded to coarser resolution by averaging all data falling into the new grid
    - Unit standardization
    - Reorder the dimensions to (lon, lat, ind)
    - Unrealistic values to NaN
    =#
    # _count = 0;
    # push!()

    return nothing
end;
=#
