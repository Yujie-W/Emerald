#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2021-Dec-24: migrate the function from PkgUtility to NetcdfIO
#     2022-Jan-28: remove the complicated funtion to create var and dim at the same time
#     2022-Jan-28: add global attributes to the generated file
#     2023-Feb-23: migrate from JuliaUtility to Emerald
#
#######################################################################################################################################################################################################
"""

    save_nc!(file::String, var_name::String, var_data::Array{T,N}, var_attribute::Dict{String,String}; compress::Int = 4, growable::Bool = false) where {T<:Union{AbstractFloat,Int,String},N}

Save the 1D, 2D, or 3D data as netcdf file, given
- `file` Path to save the dataset
- `var_name` Variable name for the data in the NC file
- `var_data` Data to save
- `var_attribute` Variable attributes for the data, such as unit and long name
- `compress` Compression level fro NetCDF, default is 4
- `growable` If true, make index growable, default is false

Note that this is a wrapper function of create_nc and append_nc:
- If var_data is 1D, the dim is set to ind
- If var_data is 2D, the dims are set to lon and lat
- If var_data is 3D, the dims are set to long, lat, and ind

#
    save_nc!(file::String, df::DataFrame, var_names::Vector{String}, var_attributes::Vector{Dict{String,String}}; compress::Int = 4, growable::Bool = false)
    save_nc!(file::String, df::DataFrame; compress::Int = 4, growable::Bool = false)

Save DataFrame to NetCDF, given
- `file` Path to save the data
- `df` DataFrame to save
- `var_names` The label of data in DataFrame to save
- `var_attributes` Variable attributes for the data to save
- `compress` Compression level fro NetCDF, default is 4
- `growable` If true, make index growable, default is false

---
# Examples
```julia
# save 1D, 2D, and 3D data
data1 = rand(12) .+ 273.15;
data2 = rand(36,18) .+ 273.15;
data3 = rand(36,18,12) .+ 273.15;

save_nc!("data1.nc", "data1", data1, Dict("description" => "Random temperature", "unit" => "K"));
save_nc!("data2.nc", "data2", data2, Dict("description" => "Random temperature", "unit" => "K"));
save_nc!("data3.nc", "data3", data3, Dict("description" => "Random temperature", "unit" => "K"));

# save DataFrame
df = DataFrame();
df[!,"A"] = rand(5);
df[!,"B"] = rand(5);
df[!,"C"] = rand(5);
save_nc!("dataf.nc", df, ["A","B"], [Dict("A" => "Attribute A"), Dict("B" => "Attribute B")]);

save_nc!("test.nc", df);
```

"""
function save_nc! end

save_nc!(file::String, var_name::String, var_data::Array{T,N}, var_attribute::Dict{String,String}; compress::Int = 4, growable::Bool = false) where {T<:Union{AbstractFloat,Int,String},N} = (
    @assert 1 <= N <= 3 "Variable must be a 1D, 2D, or 3D dataset!";
    @assert 0 <= compress <= 9 "Compression rate must be within 0 to 9";

    # create the file
    _dset = Dataset(file, "c");

    # global title attribute
    for (_title,_notes) in ATTR_ABOUT
        _dset.attrib[_title] = _notes;
    end;

    # the case if the dimension is 1D
    if N==1
        _n_ind = (growable ? Inf : length(var_data));
        _inds  = collect(eachindex(var_data));
        add_nc_dim!(_dset, "ind", _n_ind);
        append_nc!(_dset, "ind", _inds, ATTR_CYC, ["ind"]; compress=compress);
        append_nc!(_dset, var_name, var_data, var_attribute, ["ind"]; compress=compress);

        close(_dset);

        return nothing
    end;

    # if the dimension is 2D or 3D
    _n_lon   = size(var_data, 1);
    _n_lat   = size(var_data, 2);
    _res_lon = 360 / _n_lon;
    _res_lat = 180 / _n_lat;
    _lons    = collect(_res_lon/2:_res_lon:360) .- 180;
    _lats    = collect(_res_lat/2:_res_lat:180) .- 90;
    add_nc_dim!(_dset, "lon", _n_lon);
    add_nc_dim!(_dset, "lat", _n_lat);
    append_nc!(_dset, "lon", _lons, ATTR_LON, ["lon"]; compress=compress);
    append_nc!(_dset, "lat", _lats, ATTR_LAT, ["lat"]; compress=compress);

    if N==2
        append_nc!(_dset, var_name, var_data, var_attribute, ["lon", "lat"]; compress=compress);
    elseif N==3
        _n_ind = (growable ? Inf : size(var_data,3));
        _inds  = collect(1:_n_ind);
        add_nc_dim!(_dset, "ind", _n_ind);
        append_nc!(_dset, "ind", _inds, ATTR_CYC, ["ind"]; compress=compress);
        append_nc!(_dset, var_name, var_data, var_attribute, ["lon", "lat", "ind"]; compress=compress);
    end;

    close(_dset);

    return nothing
);

save_nc!(file::String, df::DataFrame, var_names::Vector{String}, var_attributes::Vector{Dict{String,String}}; compress::Int = 4, growable::Bool = false) = (
    @assert 0 <= compress <= 9 "Compression rate must be within 0 to 9";
    @assert length(var_names) == length(var_attributes) "Variable name and attributes lengths must match!";

    # create the file
    _dset = Dataset(file, "c");

    # global title attribute
    for (_title,_notes) in ATTR_ABOUT
        _dset.attrib[_title] = _notes;
    end;

    # define dimension related variables
    _n_ind = (growable ? Inf : size(df)[1]);
    _inds  = collect(1:_n_ind);

    # save the variables
    add_nc_dim!(_dset, "ind", _n_ind);
    append_nc!(_dset, "ind", _inds, ATTR_CYC, ["ind"]; compress=compress);
    for _i in eachindex(var_names)
        append_nc!(_dset, var_names[_i], df[:, var_names[_i]], var_attributes[_i], ["ind"]; compress = compress);
    end;

    close(_dset);

    return nothing
);

save_nc!(file::String, df::DataFrame; compress::Int = 4, growable::Bool = false) = (
    _var_names = names(df);
    _var_attrs = [Dict{String,String}(_vn => _vn) for _vn in _var_names];

    save_nc!(file, df, _var_names, _var_attrs; compress=compress, growable = growable);

    return nothing
);
