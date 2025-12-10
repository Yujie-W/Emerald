"""

    site_result_tuple(spac::BulkSPAC{FT}, wd::Dict{String,Vector{FT}}, saving_setting::Vector{ParameterFunctionMapper}) where {FT}

Create a NamedTuple to store simulation results, given
- `spac` the SPAC model
- `wd` the weather driver data
- `saving_setting` the dictionary to store the settings for saving the outputs

"""
function site_result_tuple(spac::BulkSPAC{FT}, wd::Dict{String,Vector{FT}}, saving_setting::Vector{ParameterFunctionMapper}) where {FT}
    # length of the time series
    n = length(wd["FDOY"]);

    # add the fields to store outputs
    new_wdf_cols = String[];
    for lpm in saving_setting
        if lpm.to_save
            push!(new_wdf_cols, lpm.name);
        end;
    end;

    return NamedTuple{Tuple(Symbol.(new_wdf_cols))}(Tuple([ones(FT,n) .* FT(NaN) for _ in new_wdf_cols]))
end;
