#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: migrate from research repo to Emerald
#     2024-Apr-15: skip regrid if the input file does not exist
#
#######################################################################################################################################################################################################
"""

    regrid_ERA5!(year::Int, nx::Int = 1; notification::Bool = false)
    regrid_ERA5!(year::Int, nx::Int, label::String, var_name::String; folder::String = ERA5_FOLDER_HOURLY)

Regrid the ERA5 datasets, given
- `year` Which year of data to regrid
- `nx` The spatial resolution is `1/nx` degree
- `notification` If true, send out emails. Default is `false`
- `label` File name label of the source NC dataset
- `var_name` Variable label in the NC dataset
- `folder` The folder where the original NC datasets are stored. Default is `ERA5_FOLDER_HOURLY`

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
function regrid_ERA5! end;

regrid_ERA5!(year::Int, nx::Int = 1; notification::Bool = false) = (
    era5_wd = ERA5SingleLevelsDriver();
    era5_labs = [getfield(era5_wd, fn)[2] for fn in fieldnames(ERA5SingleLevelsDriver)];
    era5_vars = [getfield(era5_wd, fn)[1] for fn in fieldnames(ERA5SingleLevelsDriver)];
    regrid_ERA5!.(year, nx, era5_labs, era5_vars; folder = ERA5_FOLDER_HOURLY);
    @info "Finished regridding all the datasets!";
    if notification
        send_email!("[ERA5 DATA STATUS] Regridding data for year $(year)",
                    "fluo@gps.caltech.edu",
                    "jesiner@gmail.com",
                    "ERA5 data regridding is finished for year $(year)!");
    end;

    return nothing;
);

regrid_ERA5!(year::Int, nx::Int, label::String, var_name::String; folder::String = ERA5_FOLDER_HOURLY) = (
    file_in = original_file_path(label, year; folder = folder);
    file_out = reprocessed_file_path(label, year, nx; folder = folder);

    # if file exists already, skip
    if isfile(file_out)
        @info "File $(file_out) already exists!";
        return nothing;
    end;

    # if original file does not exist, skip
    if !isfile(file_in)
        @warn "File $(file_in) does not exist, skipping...";
        return nothing;
    end;

    # read the file per slice
    @info "Reading and regridding file $(file_in) per time slice...";
    times = read_nc(file_in, "time");
    matn = zeros((360nx,180nx,length(times))) .* NaN;
    @showprogress for _itim in eachindex(times)
        mati = read_nc(file_in, var_name, _itim);
        z_ss = Int(size(mati,1) / 360);
        dxy  = Int(z_ss/nx);
        nvar = replace(mati, missing=>NaN);
        for ilon in axes(matn,1), ilat in axes(matn,2)
            # from +0 to +180
            if size(matn,1)/2 < ilon <= size(matn,1)
                ilon_raw = ilon - Int(size(matn,1)/2);
                subvar = nvar[dxy*(ilon_raw-1)+1:dxy*ilon_raw+1, dxy*(ilat-1)+1:dxy*ilat+1];
            # from -180 to -0.NNN
            elseif ilon < size(matn,1)/2
                ilon_raw = ilon + Int(size(matn,1)/2);
                subvar = nvar[dxy*(ilon_raw-1)+1:dxy*ilon_raw+1, dxy*(ilat-1)+1:dxy*ilat+1];
            # -0
            else
                ilon_raw = size(matn,1);
                subvar = [collect(nvar[dxy*(ilon_raw-1)+1:dxy*ilon_raw, dxy*(ilat-1)+1:dxy*ilat+1]); collect(nvar[1:1, dxy*(ilat-1)+1:dxy*ilat+1])];
            end;
            # for soil water content, replace with NaN if the SWC <= 0.01
            if var_name in ["SWC_1", "SWC_2", "SWC_3", "SWC_4"]
                @. subvar[subvar .<= 0] = NaN;
            end;
            matn[ilon,size(matn,2)+1-ilat,_itim] = nanmean(subvar);
        end;
    end;

    # save the regridded dataset
    @info "Saving regridded dataset to $(file_out)...";
    attr = Dict(var_name => label, "unit" => "Same as $(file_in)");
    save_nc!(file_out, var_name, matn, attr);

    # set the variables to nothing to clean the memory
    times = nothing;
    matn = nothing;
    attr = nothing;

    return nothing;
);
