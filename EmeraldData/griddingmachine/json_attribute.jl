#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-May-12: move function from GriddingMachineDatasets
#
#######################################################################################################################################################################################################
"""

    reference_attribute_dict()

Create a Dict that stores reference information about variable attributes

"""
function reference_attribute_dict()
    @info "These inputs are meant to generate the reference attributes witin the Netcdf dataset...";

    # loop the inputs until satisfied
    _attribute_dict = Dict{String,Any}();
    while true
        _msg = "    Please input the author information (e.g., Name S. et al.) > ";
        _authors = input_string(_msg; no_space = false);

        _msg = "    Please input the year of the publication > ";
        _year_pub = input_integer(_msg);

        _msg = "    Please input the title of the publication > ";
        _title = input_string(_msg; no_space = false);

        _msg = "    Please input the journal of the publication > ";
        _journal = input_string(_msg; no_space = false);

        _msg = "    Please input the DOI of the publication > ";
        _doi = input_string(_msg; no_space = false);

        _attribute_dict = Dict{String,String}(
            "authors"   => _authors,
            "year"      => _year_pub,
            "title"     => _title,
            "journal"   => _journal,
            "doi"       => _doi,
        );
        @show _attribute_dict

        # ask if the Dict looks okay, if so break
        _msg = "Is the generated dict okay? If not, type <N/n or No> to redo the inputs > ";
        _satisfy = input_yes_or_no(_msg; bool_conversion = true);
        if _satisfy
            break;
        end;
    end;

    return _attribute_dict
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-May-12: move function from GriddingMachineDatasets
#
#######################################################################################################################################################################################################
"""

    variable_attribute_dict()

Create a Dict that stores information about variable attributes

"""
function variable_attribute_dict()
    @info "These inputs are meant to generate the variable attributes witin the Netcdf dataset...";

    # loop the inputs until satisfied
    _attribute_dict = Dict{String,String}();
    while true
        _msg = "    Please input the long name of the variable to save > ";
        _longname = input_string(_msg; no_space = false);

        _msg = "    Please input the unit of the variable to save > ";
        _unit = input_string(_msg; no_space = false);

        _msg = "    Please input some more details of the variable to save > ";
        _about = input_string(_msg; no_space = false);

        _attribute_dict = Dict{String,String}(
            "long_name" => _longname,
            "unit"      => _unit,
            "about"     => _about,
        );
        @show _attribute_dict

        # ask if the Dict looks okay, if so break
        _msg = "Is the generated dict okay? If not, type <N/n or No> to redo the inputs > ";
        _satisfy = input_yes_or_no(_msg; bool_conversion = true);
        if _satisfy
            break;
        end;
    end;

    return _attribute_dict
end
