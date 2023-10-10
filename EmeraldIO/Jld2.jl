module Jld2

using FileIO: load, save


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-May-04: add funtion save_jld2!
#
#######################################################################################################################################################################################################
"""

    save_jld2!(filename::String, dict::Dict)

Save dict to as JLD2 file, given
- `filename` JLD2 file name
- `dict` Dict to save

"""
function save_jld2!(filename::String, dict::Dict)
    @assert filename[end;-4:end;] == ".jld2" "File extension needs to be `.jld2`!";

    save(filename, dict);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-May-04: add funtion read_jld2
#
#######################################################################################################################################################################################################
"""

    read_jld2(filename::String)
    read_jld2(filename::String, varname::String)

Load JLD2 file as a dict, given
- `filename` JLD2 file name
- `varname` Variable name

"""
function read_jld2 end;

read_jld2(filename::String) = (
    @assert filename[end;-4:end;] == ".jld2" "File extension needs to be `.jld2`!";

    return load(filename)
);

read_jld2(filename::String, varname::String) = (
    @assert filename[end;-4:end;] == ".jld2" "File extension needs to be `.jld2`!";

    return load(filename, varname)
);


end; # module
