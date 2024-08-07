#=
#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-May-12: move function from GriddingMachineDatasets
#
#######################################################################################################################################################################################################
"""

    variable_dict()

Return a dict of the input variable setup

"""
function variable_dict()
    # loop the inputs until satisfied
    _data_dict = Dict{String,Any}();
    while true
        _msg = "        Please input label or band name > ";
        _data_name = input_string(_msg; no_space = true);

        _msg = "        What is the longitude axis number of the data? (e.g., 1 for [lon,lat,time] and 2 for [lat,lon,time]) > ";
        _i_lon = input_integer(_msg; int_conversion = true);

        _msg = "        What is the latitude axis number of the data? (e.g., 2 for [lon,lat,time] and 1 for [lat,lon,time]) > ";
        _i_lat = input_integer(_msg; int_conversion = true);

        _msg = "        What is the index axis number of the data? (e.g., 3 for time in [lon,lat,time] and 0 for [lat,lon]) > "
        _i_ind = input_integer(_msg; int_conversion = true);
        if _i_ind == 0
            _i_ind = nothing;
        end;

        _msg = "        If you have extra scaling you want to make, type it here (NCDatasets may do that already, need to double check, example: x -> log(x)) > ";
        println(_msg);
        _scaling_function = readline();

        _msg = "        What are your mask function for NaN, type it here, e.g., x -> (0.1 < x <= 0.2 && x * 6 > 1) > ";
        println(_msg);
        _masking_function = readline();

        _data_dict = Dict{String,Any}(
            "DATA_NAME"            => _data_name,
            "LONGITUDE_AXIS_INDEX" => _i_lon,
            "LATITUDE_AXIS_INDEX"  => _i_lat,
            "INDEX_AXIS_INDEX"     => _i_ind,
            "SCALING_FUNCTION"     => _scaling_function,
            "MASKING_FUNCTION"     => _masking_function,
        );
        @show _data_dict

        # ask if the Dict looks okay, if so break
        _msg = "Is the generated dict okay? If not, type <N/n or No> to redo the inputs > ";
        _satisfy = input_yes_or_no(_msg; bool_conversion = true);
        if _satisfy
            break;
        end;
    end;

    return _data_dict
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-May-12: move function from GriddingMachineDatasets
#
#######################################################################################################################################################################################################
"""

    variable_dicts()

Create a vector of Dicts that store variable information

"""
function variable_dicts()
    @info "These questions are about how to read the data, please be careful about them...";

    # ask for how many independent variables do you want to save as DATA
    _msg = "    How many variables do you want to save as DATA (e.g., combine data1 and data2 in one netcdf file) > ";
    _data_count = input_integer(_msg; int_conversion = true);
    _data_dicts = Dict{String,Any}[];
    for _i_data in 1:_data_count
        println("    For data $(_i_data):");
        push!(_data_dicts, variable_dict());
    end;

    return _data_dicts
end;
=#
