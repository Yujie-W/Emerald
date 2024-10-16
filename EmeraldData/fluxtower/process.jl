#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-02-05: add function to process the data
#
#######################################################################################################################################################################################################
"""

    process_data!(ftds::FluxTowerDataset; displaying::Bool = false, force::Bool = false, saving::Bool = true)

Process the flux tower data and save it as a netcdf file, given
- `ftds` flux tower dataset struct
- `displaying` whether to display the processing progress
- `force` whether to force the processing (and saving)
- `saving` whether to save the processed data as a netcdf file (if not, return the processed data frame)

"""
function process_data!(ftds::FluxTowerDataset; displaying::Bool = false, force::Bool = false, saving::Bool = true)
    # file name for the processed data
    fn = "$(ftds.LABEL)_$(ftds.TRESO)_$(ftds.VER_TAG)";
    file_in = "";
    file_out = "";
    if ftds.LABEL[1:3] == "AMF"
        file_in = "$(AMERIFLUX_DATA)/$(fn).csv";
        file_out = "$(AMERIFLUX_REPROCESSED)/$(fn).nc";
    elseif ftds.LABEL[1:3] == "FLX"
        file_in = "$(FLUXNET_DATA)/$(fn).csv";
        file_out = "$(FLUXNET_REPROCESSED)/$(fn).nc";
    else
        return error("Unknown tower label: $(ftds.LABEL)")
    end;

    # if the file exists and force is not set, return nothing
    if isfile(file_out) && !force
        return nothing
    end;

    # otherwise, process the data and save the data as a netcdf file
    df_in = read_csv(file_in; skiprows = 2);
    df_out = DataFrame(
                YEAR   = [parse(Int, string(dfr.TIMESTAMP_START)[1:4]) for dfr in eachrow(df_in)],
                MONTH  = [parse(Int, string(dfr.TIMESTAMP_START)[5:6]) for dfr in eachrow(df_in)],
                DAY    = [parse(Int, string(dfr.TIMESTAMP_START)[7:8]) for dfr in eachrow(df_in)],
                HOUR   = [parse(Int, string(dfr.TIMESTAMP_START)[9:10]) for dfr in eachrow(df_in)],
                MINUTE = [parse(Int, string(dfr.TIMESTAMP_START)[11:12]) + 15 for dfr in eachrow(df_in)],
                WIND   = detect_data(df_in, "WS"; displaying = displaying),         # m s⁻¹
                VPD    = detect_data(df_in, "VPD"; displaying = displaying) .* 100, # Pa
                T_AIR  = detect_data(df_in, "TA"; displaying = displaying) .+ T₀(), # K
                RAD    = detect_data(df_in, "SW_IN"; displaying = displaying),      # W m⁻²
                RAD_LW = detect_data(df_in, "LW_IN"; displaying = displaying),      # W m⁻²
                P_ATM  = detect_data(df_in, "PA"; displaying = displaying) .* 1000, # Pa
                PRECIP = detect_data(df_in, "P"; displaying = displaying),          # m
                CO2    = detect_data(df_in, "CO2"; displaying = displaying),        # ppm
                T_SOIL = detect_data(df_in, "TS"; displaying = displaying) .+ T₀(), # K
                SWC    = detect_data(df_in, "SWC"; displaying = displaying),        # m³ m⁻³
                SW_OUT = detect_data(df_in, "SW_OUT"; displaying = displaying),     # W m⁻²
                LW_OUT = detect_data(df_in, "LW_OUT"; displaying = displaying),     # W m⁻²
                NEE    = detect_data(df_in, "FC"; displaying = displaying) .+
                         detect_data(df_in, "SC"; displaying = displaying),         # µmol m⁻² s⁻¹
    );

    # if saving is set, save the data as a netcdf file
    if saving
        @info "Saving the processed data to $(file_out)";
        save_nc!(file_out, df_out);
        return nothing
    else
        return df_out
    end;
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-02-05: add function to detect the data from input data frame
#
#######################################################################################################################################################################################################
"""

    detect_data(df::DataFrame, data_label::String; displaying::Bool = false)

Detect the data from the input data frame, given
- `df` input data frame
- `data_label` label of the column that start with the `data_label` string
- `displaying` whether to display the processing progress

"""
function detect_data(df::DataFrame, data_label::String; displaying::Bool = false)
    label = findall(x -> startswith(x, "$(data_label)_") && !startswith(x, "$(data_label)_PI_"), names(df));
    label_pi = findall(x -> startswith(x, "$(data_label)_PI_") && !startswith(x, "$(data_label)_PI_F_"), names(df));
    label_pi_f = findall(x -> startswith(x, "$(data_label)_PI_F_"), names(df));
    data = df[:, label];
    data_pi = df[:, label_pi];
    data_pi_f = df[:, label_pi_f];

    if displaying
        @info "The labels associated with $(data_label) are:";
        @show names(df)[label] names(df)[label_pi] names(df)[label_pi_f];
    end;

    # loop through the data and find the first non-NaN mean value
    output = Float64[];
    for i in axes(df,1)
        vec_data = length(label) >= 1 ? collect(Float64, data[i,:]) : Float64[];
        vec_data_pi = length(label_pi) >= 1 ? collect(Float64, data_pi[i,:]) : Float64[];
        vec_data_pi_f = length(label_pi_f) >= 1 ? collect(Float64, data_pi_f[i,:]) : Float64[];

        # if the data == -9999, replace it with NaN
        @. vec_data[vec_data .<= -9990] = NaN;
        @. vec_data_pi[vec_data_pi .<= -9990] = NaN;
        @. vec_data_pi_f[vec_data_pi_f .<= -9990] = NaN;

        # take the nanmean of each vector
        mean_data = nanmean(vec_data);
        mean_data_pi = nanmean(vec_data_pi);
        mean_data_pi_f = nanmean(vec_data_pi_f);

        # if mean_data is not NaN, use it as the data
        if !isnan(mean_data)
            push!(output, mean_data);
        # if mean_data is NaN but mean_data_pi is not NaN, use it as the data
        elseif !isnan(mean_data_pi)
            push!(output, mean_data_pi);
        # if mean_data is NaN but mean_data_pi is NaN but mean_data_pi_f is not NaN, use it as the data
        elseif !isnan(mean_data_pi_f)
            push!(output, mean_data_pi_f);
        # otherwise, use NaN as the data
        else
            push!(output, NaN);
        end;
    end;

    return output
end;
