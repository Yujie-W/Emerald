module EmeraldEarth

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using JLD2: save

using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: lat_ind, lon_ind, read_LUT

using ..EmeraldIO.Text: read_csv
using ..EmeraldMath.Stats: nanmax
using ..EmeraldCore

# simulation settings
RESULT_FOLDER = "/home/wyujie/DATASERVER/model/CLIMA/LAND/simulations";
SETUP_FOLDER  = "/home/wyujie/DATASERVER/model/CLIMA/LAND/setups";

# CLM5 settings
CLM5_PFTG = [0, 2.35, 2.35, 2.35, 4.12, 4.12, 4.45, 4.45, 4.45, 4.7, 4.7, 4.7, 2.22, 5.25, 1.62, 5.79, 5.79] .* sqrt(1000);
CLM5_PFTS = ["not_vegetated",
             "needleleaf_evergreen_temperate",
             "needleleaf_evergreen_boreal",
             "needleleaf_deciduous_boreal",
             "broadleaf_evergreen_tropical",
             "broadleaf_evergreen_temperate",
             "broadleaf_deciduous_tropical",
             "broadleaf_deciduous_temperate",
             "broadleaf_deciduous_boreal",
             "evergreen_shrub",
             "deciduous_temperate_shrub",
             "deciduous_boreal_shrub",
             "c3_arctic_grass",
             "c3_non-arctic_grass",
             "c4_grass",
             "c3_crop",
             "c3_irrigated"];


include("griddingmachine.jl")


end # module
