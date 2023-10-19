#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-May-12: move function from GriddingMachineDatasets
#
#######################################################################################################################################################################################################
"""

    griddingmachine_dict()

Create a Dict that stores information about GriddingMachine tag

"""
function griddingmachine_dict()
    @info "These inputs are meant to generate the GriddingMachine TAG...";

    # loop the inputs until satisfied
    _griddingmachine_dict = Dict{String,Any}();
    while true
        _msg = "    Please indicate the level 1 label for the dataset (e.g., GPP as in GPP_VPM_2X_1M_V1) > ";
        _label = input_string(_msg, uppercase; no_space = true);

        _msg = "    Please indicate the level 2 label for the dataset (e.g., VPM as in GPP_VPM_2X_1M_V1, leave empty is there is not any) > ";
        _label_extra = input_string(_msg, to_nothing_or_uppercase; allow_blank = true, no_space = true);

        _msg = "    Please indicate the spatial resolution represented with an integer (N for 1/N °) > ";
        _spatial_resolution_nx = input_integer(_msg; int_conversion = true, non_negative = true);

        _msg = "    Please indicate the temporal resolution (e.g., 8D, 1M, and 1Y) > ";
        _temporal_resolution = verified_input(_msg, uppercase, is_gm_mt);

        _msg = "    Please indicate the range of years (e.g., 2001:2022, and 2001,2005,2020 and empty for non-specific) > ";
        _years = verified_input(_msg, to_nothing_or_int_vector, is_nothing_or_vector);

        _msg = "    Please indicate the version number of the dataset (1 for V1) > ";
        _version = input_integer(_msg; int_conversion = true, positive = true);

        _labeling = isnothing(_label_extra) ? _label : _label * "_" * _label_extra;
        if isnothing(_years)
            @info "The GriddingMachine tag will be $(_labeling)_$(_spatial_resolution_nx)X_$(_temporal_resolution)_V$(_version)";
        else
            @info "The GriddingMachine tag will be $(_labeling)_$(_spatial_resolution_nx)X_$(_temporal_resolution)_YEAR_V$(_version)";
        end;

        # ask if the TAG looks okay, if so break
        _msg = "Is the generated tag okay? If not, type <N/n or No> to redo the inputs > ";
        _satisfy = input_yes_or_no(_msg; bool_conversion = true);
        if _satisfy
            _griddingmachine_dict = Dict{String,Any}(
                "LABEL"         => _label,
                "EXTRA_LABEL"   => _label_extra,
                "LAT_LON_RESO"  => _spatial_resolution_nx,
                "TEMPORAL_RESO" => _temporal_resolution,
                "YEARS"         => _years,
                "VERSION"       => _version,
            );
            break;
        end;
    end;

    return _griddingmachine_dict
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-May-12: move function from GriddingMachineDatasets
#
#######################################################################################################################################################################################################
"""

    griddingmachine_tag(dict::Dict, year::Int = 0)

Generate the GriddingMachine TAG, given
- `dict` Dict readed from JSON file
- `year` Which year (used only when YEARS in dict is not nothing)

"""
function griddingmachine_tag(dict::Dict, year::Int = 0)
    _griddingmachine_dict = dict["GRIDDINGMACHINE"];

    _years = _griddingmachine_dict["YEARS"];
    _label = _griddingmachine_dict["LABEL"];
    _label_extra = _griddingmachine_dict["EXTRA_LABEL"];
    _labeling = isnothing(_label_extra) ? _label : _label * "_" * _label_extra;
    _spatial_resolution_nx = _griddingmachine_dict["LAT_LON_RESO"];
    _temporal_resolution = _griddingmachine_dict["TEMPORAL_RESO"];
    _version = _griddingmachine_dict["VERSION"];

    if isnothing(_years)
        _tag = "$(_labeling)_$(_spatial_resolution_nx)X_$(_temporal_resolution)_V$(_version)";
    else
        _tag = "$(_labeling)_$(_spatial_resolution_nx)X_$(_temporal_resolution)_$(year)_V$(_version)";
    end;

    return _tag
end;
