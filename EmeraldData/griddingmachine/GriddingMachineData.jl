#=
module GriddingMachineData

using JSON;

using NetcdfIO: append_nc!, read_nc, save_nc!

#using ..EmeraldIO.Terminal:
using ..EmeraldIO.Terminal: has_no_space, input_integer, input_string, input_yes_or_no, verified_input
using ..EmeraldUtility.Artifact: deploy_artifact!


GRIDDING_MACHINE_HOME = joinpath(homedir(), "GriddingMachine");
ARTIFACT_TOML         = "$(GRIDDING_MACHINE_HOME)/Artifacts.toml";
DATASET_FOLDER        = "$(GRIDDING_MACHINE_HOME)/reprocessed";
ARTIFACT_FOLDER       = "$(GRIDDING_MACHINE_HOME)/artifacts"
FTP_URLS              = ["ftp://fluo.gps.caltech.edu/XYZT_GRIDDING_MACHINE/artifacts"];


include("data_read.jl");
include("data_reprocess.jl");
include("deploy.jl");
include("json_attribute.jl");
include("json_data.jl");
include("json_griddingmachine.jl");
include("json_map.jl");
include("json_save.jl");
include("terminal.jl");


end; # module
=#
