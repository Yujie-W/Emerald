module GriddingMachineData

using JSON;

using ..EmeraldIO.Netcdf: append_nc!, read_nc, save_nc!
using ..EmeraldIO.Terminal: verified_input
using ..EmeraldUtility.Artifact: deploy_artifact!
using ..EmeraldVisualization: animate_data!


GRIDDING_MACHINE_HOME = "/home/wyujie/GriddingMachine";
ARTIFACT_TOML         = "$(GRIDDING_MACHINE_HOME)/Artifacts.toml";
DATASET_FOLDER        = "$(GRIDDING_MACHINE_HOME)/reprocessed";
ARTIFACT_FOLDER       = "$(GRIDDING_MACHINE_HOME)/artifacts"
FTP_URLS              = ["ftp://fluo.gps.caltech.edu/XYZT_GRIDDING_MACHINE/artifacts"];


include("griddingmachine/data_read.jl");
include("griddingmachine/data_reprocess.jl");
include("griddingmachine/deploy.jl");
include("griddingmachine/json_attribute.jl");
include("griddingmachine/json_data.jl");
include("griddingmachine/json_griddingmachine.jl");
include("griddingmachine/json_map.jl");
include("griddingmachine/json_save.jl");


end # module
