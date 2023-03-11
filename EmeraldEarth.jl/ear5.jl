#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: migrate from research repo to Emerald
#
#######################################################################################################################################################################################################
"""

    fetch_ERA5_data!(year::Int; notification::Bool = false)

Fetch all the ERA5 data rquired to run this project, given
- `year` Which year of data to download
- `notification` If true, send out emails. Default is `false`

Note that ERA5 single levels data ranges from 1979 to present, and ERA5 land data ranges from 1981 to present. Be aware not to download data out of the range.

"""
function fetch_ERA5_data!(year::Int; notification::Bool = false)
    _dts = ERA5SingleLevelsHourly();
    _dir = "$(ERA5_FOLDER)/original/";
    if !isdir(_dir)
        @warn "$(_dir) does not exist, skipping...";
        return nothing
    end

    # An email will be sent out per year, comment it if you do not want
    fetch_data!(_dts, year; vars=ERA5_LABELS, folder=_dir);
    @info "Finished downloading all the datasets!";
    if notification
        #=
        send_email!("[ERA5 DATA STATUS] Downloading data for year $(year)",
                    "fluo@gps.caltech.edu",
                    "jesiner@gmail.com",
                    "ERA5 data downloading is finished for year $(year)!");
        =#
    end;

    return nothing
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: migrate from research repo to Emerald
#
#######################################################################################################################################################################################################
"""

    regrid_ERA5!(year::Int, zoom::Int = 1; notification::Bool = false)
    regrid_ERA5!(year::Int, zoom::Int, label::String, var_name::String)

Regrid the ERA5 datasets, given
- `year` Which year of data to regrid
- `zoom` The spatial resolution is `1/zoom` degree
- `notification` If true, send out emails. Default is `false`
- `label` File name label of the source NC dataset
- `var_name` Variable label in the NC dataset

To reduce memory allocation, `regrid_ERA5!` reads in the data per slice (in time index) and regrid the data per slice. What `regrid_ERA5!` does are

- determine if the output file exists
- if `true`, skip and do nothing
- if `false`
  - read the data per slice
  - replace missing values with NaNs
  - use mean of the target region to fill output matrix
  - save the regridded matrix to a new NC file

Note that ERA5 NC datasets differ from `GriddingMachine.jl` standard in that

- Latitude is from 90 to -90 in ERA5 NC dataset, and has `180N+1` elements
- Longitude is from 0 to 359.XXX in ERA5 NC dataset, and has `360N` elements

Thus, we need to regrid the dataset to ensure that the pixel orders match. For example, if the source dataset spatial resolution is `4X` and the target spatial resolution is `1X`, we need take the
    average of data in the pixel of 120E-121E (the steps are 120, 120.25, 120.5, 120.75, 121), 30N-31N (the steps are 30, 30.25, 30.5, 30.75, 31), namely 25 elements in total. The special case is
    when the longitude ranges from 179E to 180E, and in this case, we need to include the -180E slice. Otherwise, the mean longitude is 179.375 rather than 179.5. The general function is

"""
function regrid_ERA5! end

regrid_ERA5!(year::Int, zoom::Int = 1; notification::Bool = false) = (
    regrid_ERA5!.(year, zoom, ERA5_LABELS, ERA5_LAYERS);
    @info "Finished regridding all the datasets!";
    if notification
        #=
        send_email!("[ERA5 DATA STATUS] Regridding data for year $(year)",
                    "fluo@gps.caltech.edu",
                    "jesiner@gmail.com",
                    "ERA5 data regridding is finished for year $(year)!");
        =#
    end;

    return nothing;
);

regrid_ERA5!(year::Int, zoom::Int, label::String, var_name::String) = (
    _file_in  = "$(ERA5_FOLDER)/original/$(label)_SL_$(year).nc";
    _file_out = "$(ERA5_FOLDER)/reprocessed/$(label)_SL_$(year)_$(zoom)X.nc";

    # if file exists already, skip
    if isfile(_file_out)
        @info "File $(_file_out) already exists!";
        return nothing;
    end;

    # read the file per slice
    @info "Reading and regridding file $(_file_in) per time slice...";
    _time = read_nc(_file_in, "time");
    _matn = zeros((360*zoom,180*zoom,length(_time))) .* NaN;
    @showprogress for _itim in eachindex(_time)
        _mati = read_nc(_file_in, var_name, _itim);
        _z_ss = Int(size(_mati,1) / 360);
        _dxy  = Int(_z_ss/zoom);
        _nvar = replace(_mati, missing=>NaN);
        for _ilon in axes(_matn,1), _ilat in axes(_matn,2)
            # from +0 to +180
            if size(_matn,1)/2 < _ilon <= size(_matn,1)
                __ilon = _ilon - Int(size(_matn,1)/2);
                _sub   = _nvar[_dxy*(__ilon-1)+1:_dxy*__ilon+1, _dxy*(_ilat-1)+1:_dxy*_ilat+1];
            # from -180 to -0.NNN
            elseif _ilon < size(_matn,1)/2
                __ilon = _ilon + Int(size(_matn,1)/2);
                _sub   = _nvar[_dxy*(__ilon-1)+1:_dxy*__ilon+1, _dxy*(_ilat-1)+1:_dxy*_ilat+1];
            # -0
            else
                __ilon = size(_matn,1);
                _sub   = [collect(_nvar[_dxy*(__ilon-1)+1:_dxy*__ilon, _dxy*(_ilat-1)+1:_dxy*_ilat+1]); collect(_nvar[1:1, _dxy*(_ilat-1)+1:_dxy*_ilat+1])];
            end;
            # for soil water content, replace with NaN if the SWC <= 0.01
            if var_name in ["SWC_1", "SWC_2", "SWC_3", "SWC_4"]
                _sub[_sub .<= 0.01] .= NaN;
            end;
            _matn[_ilon,size(_matn,2)+1-_ilat,_itim] = nanmean(_sub);
        end;
    end;

    # save the regridded dataset
    @info "Saving regridded dataset to $(_file_out)...";
    _attr = Dict(var_name => label * "_SL", "unit" => "Same as $(_file_in)");
    save_nc!(_file_out, var_name, _matn, _attr);

    return nothing;
);
