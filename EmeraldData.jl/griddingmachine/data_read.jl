"""

    read_data_2d(data::Array, ind::Int, dict::Dict, flipping::Vector; scaling_function::Union{Function,Nothing} = nothing)

Return formatted 2D data, given
- `data` Input 3D data
- `ind` Index for 3rd dimension
- `dict` Dict about the data format
- `flipping` Whether to flip latitude and longitude
- `scaling_function` Scaling function that to apply

"""
function read_data_2d(data::Array, ind::Int, dict::Dict, flipping::Vector; scaling_function::Union{Function,Nothing} = nothing)
    # read the layer based on the index orders
    if isnothing(dict["INDEX_AXIS_INDEX"])
        if dict["LONGITUDE_AXIS_INDEX"] == 1 && dict["LATITUDE_AXIS_INDEX"] == 2
            _eata = data;
        else
            _eata = data';
        end;
    else
        if dict["INDEX_AXIS_INDEX"] == 3
            if dict["LONGITUDE_AXIS_INDEX"] == 1 && dict["LATITUDE_AXIS_INDEX"] == 2
                _eata = data[:,:,ind];
            else
                _eata = data[:,:,ind]';
            end;
        elseif dict["INDEX_AXIS_INDEX"] == 2
            if dict["LONGITUDE_AXIS_INDEX"] == 1 && dict["LATITUDE_AXIS_INDEX"] == 3
                _eata = data[:,ind,:];
            else
                _eata = data[:,ind,:]';
            end;
        elseif dict["INDEX_AXIS_INDEX"] == 1
            if dict["LONGITUDE_AXIS_INDEX"] == 2 && dict["LATITUDE_AXIS_INDEX"] == 3
                _eata = data[ind,:,:];
            else
                _eata = data[ind,:,:]';
            end;
        end;
    end;

    # flip the lat and lons
    _fata = flipping[1] ? _eata[:,end:-1:1] : _eata;
    _gata = flipping[2] ? _fata[end:-1:1,:] : _fata;

    # add a scaling function
    _hata = isnothing(scaling_function) ? _gata : scaling_function.(_gata);

    return _hata
end


"""

    read_data(filename::String, dict::Dict, flipping::Vector; scaling_function::Union{Function,Nothing} = nothing)

Return the formatted data, given
- `filename` File to read
- `dict` Dict about the data format
- `flipping` Whether to flip latitude and longitude
- `scaling_function` Scaling function that to apply

"""
function read_data(filename::String, dict::Dict, flipping::Vector; scaling_function::Union{Function,Nothing} = nothing)
    _data = read_nc(filename, dict["DATA_NAME"]);

    # rotate the data if necessary
    if isnothing(dict["INDEX_AXIS_INDEX"])
        return read_data_2d(_data, 1, dict, flipping; scaling_function = scaling_function)
    else
        _eata = zeros(Float64, size(_data, dict["LONGITUDE_AXIS_INDEX"]), size(_data, dict["LATITUDE_AXIS_INDEX"]), size(_data, dict["INDEX_AXIS_INDEX"]));
        for _ind in axes(_data, dict["INDEX_AXIS_INDEX"])
            _eata[:,:,_ind] .= read_data_2d(_data, _ind, dict, flipping; scaling_function = scaling_function);
        end;

        return _eata
    end;
end
