#######################################################################################################################################################################################################
#
# Changes to the function
# General:
#     2021-Dec-24: move the function from PkgUtility to NetcdfIO
#     2022-Jan-28: remove the complicated funtion to create var and dim at the same time
#     2022-Jan-28: add method to add data to Dataset
#     2023-Feb-23: migrate from JuliaUtility to Emerald
#
#######################################################################################################################################################################################################
"""

    append_nc!(ds::Dataset, var_name::String, var_data::Array{T,N}, var_attributes::Dict{String,String}, dim_names::Vector{String}; compress::Int = 4) where {T<:Union{AbstractFloat,Int,String},N}
    append_nc!(file::String, var_name::String, var_data::Array{T,N}, var_attributes::Dict{String,String}, dim_names::Vector{String}; compress::Int = 4) where {T<:Union{AbstractFloat,Int,String},N}

Append data to existing netcdf dataset, given
- `ds` A `NCDatasets.Dataset` type dataset
- `var_name` New variable name to write to
- `var_data` New variable data to write, can be integer, float, and string with N dimens
- `var_attributes` New variable attributes
- `dim_names` Dimension names in the netcdf file
- `compress` Compression level fro NetCDF, default is 4
- `file` Path of the netcdf dataset

---
# Examples
```julia
create_nc!("test.nc", String["lon", "lat", "ind"], [36, 18, 5]);

dset = Dataset("test.nc", "a");
append_nc!(dset, "str", ["A" for i in 1:18], Dict("longname" => "test strings"), ["lat"]);
append_nc!(dset, "d3d", rand(36,18,5), Dict("longname" => "a 3d dataset"), ["lon", "lat", "ind"]);
close(dset);

append_nc!("test.nc", "lat", collect(1:18), Dict("longname" => "latitude"), ["lat"]);
append_nc!("test.nc", "d2d", rand(36,18), Dict("longname" => "a 2d dataset"), ["lon", "lat"]);
```

"""
function append_nc! end

append_nc!(ds::Dataset, var_name::String, var_data::Array{T,N}, var_attributes::Dict{String,String}, dim_names::Vector{String}; compress::Int = 4) where {T<:Union{AbstractFloat,Int,String},N} = (
    # only if variable does not exist create the variable
    @assert !(var_name in keys(ds)) "You can only add new variable to the dataset!";
    @assert length(dim_names) ==  N "Dimension must be match!";
    @assert 0 <= compress <= 9 "Compression rate must be within 0 to 9";

    # if type of variable is string, set deflatelevel to 0
    if T == String
        compress = 0;
    end;

    _ds_var = defVar(ds, var_name, T, dim_names; attrib = var_attributes, deflatelevel = compress);
    _ds_var[:,:] = var_data;

    return nothing
);

append_nc!(file::String, var_name::String, var_data::Array{T,N}, var_attributes::Dict{String,String}, dim_names::Vector{String}; compress::Int = 4) where {T<:Union{AbstractFloat,Int,String},N} = (
    _dset = Dataset(file, "a");
    append_nc!(_dset, var_name, var_data, var_attributes, dim_names; compress = compress);
    close(_dset);

    return nothing
);
