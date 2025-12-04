
function site_result_tuple(spac::BulkSPAC{FT}, wd::Dict{String,Vector{FT}}, saving_dict::Dict{String,Bool}) where {FT}
    # length of the time series
    n = length(wd["FDOY"]);

    # add the fields to store outputs
    new_wdf_cols = String[];
    if saving_dict["MOD_SWC"]
        for i in eachindex(spac.soils)
            push!(new_wdf_cols, "MOD_SWC_$i");
            push!(new_wdf_cols, "MOD_SWC_ICE_$i");
        end;
    end;
    if saving_dict["MOD_P_SOIL"]
        for i in eachindex(spac.soils)
            push!(new_wdf_cols, "MOD_P_SOIL_$i");
        end;
    end;
    if saving_dict["MOD_T_SOIL"]
        for i in eachindex(spac.soils)
            push!(new_wdf_cols, "MOD_T_SOIL_$i");
        end;
    end;
    if saving_dict["MOD_T_LEAF"]
        for i in eachindex(spac.plant.leaves)
            push!(new_wdf_cols, "MOD_T_LEAF_$i");
        end;
    end;
    if saving_dict["MOD_T_MMM"]
        push!(new_wdf_cols, "MOD_T_L_MAX");
        push!(new_wdf_cols, "MOD_T_L_MEAN");
        push!(new_wdf_cols, "MOD_T_L_MIN");
    end;
    if saving_dict["MOD_P_LEAF"]
        for i in eachindex(spac.plant.leaves)
            push!(new_wdf_cols, "MOD_P_LEAF_$i");
        end;
    end;
    if saving_dict["MOD_P_MMM"]
        push!(new_wdf_cols, "MOD_P_L_MAX");
        push!(new_wdf_cols, "MOD_P_L_MEAN");
        push!(new_wdf_cols, "MOD_P_L_MIN");
    end;
    if saving_dict["MOD_ΦDΦN"]
        push!(new_wdf_cols, "ΦD");
        push!(new_wdf_cols, "ΦN");
    end;
    if saving_dict["MOD_ΦFΦP"]
        push!(new_wdf_cols, "ΦF");
        push!(new_wdf_cols, "ΦP");
    end;
    if saving_dict["MOD_HEAT"]
        push!(new_wdf_cols, "MOD_LATENT_HEAT");
        push!(new_wdf_cols, "MOD_SENSIBLE_HEAT");
        push!(new_wdf_cols, "MOD_NET_LONGWAVE");
        push!(new_wdf_cols, "MOD_NET_SHORTWAVE");
        push!(new_wdf_cols, "MOD_LONGWAVE_OUT");
        push!(new_wdf_cols, "MOD_SHORTWAVE_OUT");
    end;
    # if the label does not contain MOD_ prefix
    for label in keys(saving_dict)
        if !occursin("MOD_", label)
            if saving_dict[label]
                push!(new_wdf_cols, label);
            end;
        end;
    end;

    return NamedTuple{Tuple(Symbol.(new_wdf_cols))}(Tuple([ones(FT,n) .* FT(NaN) for _ in new_wdf_cols]))
end;
