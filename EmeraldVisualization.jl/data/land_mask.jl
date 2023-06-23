#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Jun-23: add new function to tune data based on land mask to filter out NaNs
#
#######################################################################################################################################################################################################
"""

    land_mask_data!(data::Array{FT,N}; division::Int = 1, replacement::Number = 0, threshold::Number = 0.01) where {FT,N}

Filter out NaNs in data based on land mask, given
- `data` Data to filter out NaNs
- `division` Spatial resolution of the land mask is 1/division degree
- `replacement` Value to replace NaNs
- `threshold` Threshold of land mask to filter out NaNs

"""
function land_mask_data!(data::Array{FT,N}; division::Int = 1, replacement::Number = 0, threshold::Number = 0.01) where {FT,N}
    # read land mask data from GriddingMachine if the data is not yet loaded into memory
    if isnothing(LAND_MASK)
        global LAND_MASK = read_LUT(query_collection("LM_4X_1Y_V1"))[1];
    end;

    # regrid land mask data to the same resolution as data
    _land_mask = regrid(LAND_MASK, division);

    # filter out NaNs in data based on land mask threshold
    _nan_mask = isnan.(data) .&& (_land_mask .>= threshold);

    # replace NaNs in data with replacement
    data[_nan_mask] .= replacement;

    return nothing
end
