"""
Configuration information of the SPAC model
"""
Base.@kwdef mutable struct SPACConfigInfo
    "Data name of the JLD2 file"
    DATASET::String = OLD_PHI_2021
    "JLD2 file name"
    JLD2_FILE::String = LAND_ARTIFACT
    "Message level (0 for no, 1 for progress bar, and 2 for ind)"
    MESSAGE_LEVEL::Int = 0
end;
