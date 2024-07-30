# The design of a refactored struct consists of several major components:
#     state variables are all store in a single struct as a field
#     auxiliary variables are all store in a single struct as a field
#     other sublevel structs that has their own state and auxiliary variables are stored as fields
# To avoid any confusion, it is recommended to use only state and auxiliary variables in the main struct, or use only the fields of sublevel structs

module Namespace

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS

using ..EmeraldIO.Jld2: read_jld2
using ..EmeraldIO.Text: read_csv
using ..EmeraldMath.Solver: NewtonBisectionMethod, SolutionTolerance, find_zero
using ..EmeraldPhysics.Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, GRAVITY, M_H₂O, P_ATM, T₀, T₂₅, ρ_H₂O
using ..EmeraldUtility.StructEqual: sync_struct!


# Please do not use V1/V2/V3 files here as they do not contain the Phi_PSI and Phi_PSII variables
const LAND_ARTIFACT    = artifact"land_model_spectrum_V8" * "/land_model_spectrum_V8.jld2";
const OLD_PHI_2017     = "oldphi_2017";
const OLD_PHI_2021     = "oldphi_2021";
const NEW_PHI_2017     = "newphi_2017";
const NEW_PHI_2021     = "newphi_2021";
const OLD_PHI_2017_1NM = "oldphi_2017_1nm";
const OLD_PHI_2021_1NM = "oldphi_2021_1nm";
const NEW_PHI_2017_1NM = "newphi_2017_1nm";
const NEW_PHI_2021_1NM = "newphi_2021_1nm";
const SOIL_TEXT        = read_csv("$(@__DIR__)/../../data/SOIL-TEXTURE.csv");


# General instructions to run SPAC (dependent on config)
include("general.jl");


# General methods (for users to choose from)
include("method/beta.jl");
include("method/colimit.jl");
include("method/fluorescence.jl");
include("method/kinetics.jl");
include("method/lidf.jl");
include("method/pv.jl");
include("method/soil.jl");
include("method/soil_albedo.jl");
include("method/stomata.jl");
include("method/xylem.jl");


# The configuration of the SPAC system
include("config/spectra.jl");
include("config/trace.jl");

include("config/config.jl");


# Plant hydraulics (dependent on config and method)
include("xylem/energy.jl");

include("xylem/junction.jl");
include("xylem/xylem.jl");


# Soil
include("soil/bulk.jl");
include("soil/layer.jl");


# Root system (dependent on xylem)
include("root/rhizosphere.jl");

include("root/root.jl");


# Stem system (dependent on xylem)
include("stem/stem.jl");


# Leaf system (dependent on xylem)
include("leaf/biophysics.jl");
include("leaf/energy.jl");
include("leaf/extraxylem.jl");
include("leaf/flux.jl");
include("leaf/layerflux.jl");

include("leaf/photosynthesis.jl");

include("leaf/layer.jl");
include("leaf/leaf.jl");


# Canopy
include("canopy/sensor_geometry.jl");
include("canopy/structure.jl");
include("canopy/sun_geometry.jl");

include("canopy/canopy.jl");


# Environment
include("environment/air.jl");
include("environment/radiation.jl");

include("environment/meteorology.jl");


# SPAC
include("spac/cache.jl");
include("spac/info.jl");
include("spac/memory.jl");

include("spac/plant.jl");

include("spac/bulk.jl");


end;
