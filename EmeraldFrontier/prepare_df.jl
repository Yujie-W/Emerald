#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Aug-06: isolate the function to prepare the DataFrame for the WDF
#     2025-Jun-03: add function chunk to save Φ_D and Φ_N
#
#######################################################################################################################################################################################################
"""

    prepare_wdf(spac::BulkSPAC{FT}, df::DataFrame; saving_dict::Dict{String,Any} = SAVING_DICT) where {FT}

Prepare the DataFrame for the WDF by adding the fields to store the outputs, given
- `spac::BulkSPAC{FT}`: the SPAC model
- `df::DataFrame`: the DataFrame to store the outputs
- `saving_dict::Dict{String,Any}`: the dictionary to store the settings for saving the outputs

"""
function prepare_wdf(spac::BulkSPAC{FT}, df::DataFrame; saving_dict::Dict{String,Any} = SAVING_DICT) where {FT}
    # add the fields to store outputs
    new_df_cols = String[];
    if saving_dict["MOD_SWC"]
        for i in eachindex(spac.soils)
            push!(new_df_cols, "MOD_SWC_$i");
        end;
    end;
    if saving_dict["MOD_P_SOIL"]
        for i in eachindex(spac.soils)
            push!(new_df_cols, "MOD_P_SOIL_$i");
        end;
    end;
    if saving_dict["MOD_T_SOIL"]
        for i in eachindex(spac.soils)
            push!(new_df_cols, "MOD_T_SOIL_$i");
        end;
    end;
    if saving_dict["MOD_T_LEAF"]
        for i in eachindex(spac.plant.leaves)
            push!(new_df_cols, "MOD_T_LEAF_$i");
        end;
    end;
    if saving_dict["MOD_T_MMM"]
        push!(new_df_cols, "MOD_T_L_MAX");
        push!(new_df_cols, "MOD_T_L_MEAN");
        push!(new_df_cols, "MOD_T_L_MIN");
    end;
    if saving_dict["MOD_P_LEAF"]
        for i in eachindex(spac.plant.leaves)
            push!(new_df_cols, "MOD_P_LEAF_$i");
        end;
    end;
    if saving_dict["MOD_P_MMM"]
        push!(new_df_cols, "MOD_P_L_MAX");
        push!(new_df_cols, "MOD_P_L_MEAN");
        push!(new_df_cols, "MOD_P_L_MIN");
    end;
    if saving_dict["MOD_ΦDΦN"]
        push!(new_df_cols, "ΦD");
        push!(new_df_cols, "ΦN");
    end;
    if saving_dict["MOD_ΦFΦP"]
        push!(new_df_cols, "ΦF");
        push!(new_df_cols, "ΦP");
    end;
    # if the label does not contain MOD_ prefix
    for label in keys(saving_dict)
        if !occursin("MOD_", label)
            if saving_dict[label]
                push!(new_df_cols, label);
            end;
        end;
    end;
    for label in new_df_cols
        df[!,label] .= NaN;
    end;

    # convert the DataFrame to NamedTuple
    return NamedTuple{Tuple(Symbol.(names(df)))}(Tuple([df[:,n] for n in names(df)]));
end;
