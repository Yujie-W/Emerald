#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-02-05: add function to query the processed data
#
#######################################################################################################################################################################################################
"""

    processed_flux_tower_data(ftds::FluxTowerDataset)

Query the processed flux tower data (generate if data does not exist), given
- `ftds` a FluxTowerDataset struct

"""
function processed_flux_tower_data(ftds::FluxTowerDataset)
    # file name for the processed data
    fn = "$(ftds.LABEL)_$(ftds.TRESO)_$(ftds.VER_TAG)";
    file_out = "";
    if ftds.LABEL[1:3] == "AMF"
        file_out = "$(AMERIFLUX_REPROCESSED)/$(fn).nc";
    elseif ftds.LABEL[1:3] == "FLX"
        file_out = "$(FLUXNET2015_REPROCESSED)/$(fn).nc";
    end;

    # if the file does not exist, process the data
    if !isfile(file_out)
        process_data!(ftds; displaying = true, force = false, saving = true);
    end;

    # return the processed data as a data frame
    return read_nc(file_out)
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-02-05: add function to query the raw flux tower data
#
#######################################################################################################################################################################################################
"""

    raw_flux_tower_data(ftds::FluxTowerDataset)

Query the raw flux tower data, given
- `ftds` a FluxTowerDataset struct

"""
function raw_flux_tower_data(ftds::FluxTowerDataset)
    # file name for the processed data
    fn = "$(ftds.LABEL)_$(ftds.TRESO)_$(ftds.VER_TAG)";
    file_in = "";
    if ftds.LABEL[1:3] == "AMF"
        file_in = "$(AMERIFLUX_DATA)/$(fn).csv";
    elseif ftds.LABEL[1:3] == "FLX"
        file_in = "$(FLUXNET2015_DATA)/$(fn).csv";
    end;

    # if the file does not exist, warn the user
    if !isfile(file_in)
        return error("The file $(file_in) does not exist, please double check you FluxTowerDataset struct!")
    end;

    # return the raw data as a data frame
    return read_csv(file_in; skiprows = 2)
end;
