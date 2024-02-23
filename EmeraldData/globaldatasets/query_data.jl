#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to determine the day bounds of input GriddingMachine drivers
#
#######################################################################################################################################################################################################
"""

    griddingmachine_data(data::Union{Number,Vector}, year::Int, d::Int)

Return the index of data, given
- `data` Time series of input data (or a single value)
- `year` Year
- `d` Day number

"""
function query_griddingmachine_data(data::Union{Number,Vector}, year::Int, d::Int)
    bounds = [0, 367];
    n = length(data);
    if n == 1
        bounds = [0, 367]
    elseif n == 12
        bounds = isleapyear(year) ? MDAYS_LEAP : MDAYS;
    elseif n == 46
        bounds = [collect(0:8:361); 367]
    elseif n == 52
        bounds = [collect(0:7:361); 367]
    else
        error("This temporal resolution is not supported: $(n)!");
    end;

    ind = findfirst(d .<= bounds) - 1;

    return data[ind]
end;
