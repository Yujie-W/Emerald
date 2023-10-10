# The design of a refactored struct consists of several major components:
#     state variables are all store in a single struct as a field
#     auxiliary variables are all store in a single struct as a field
#     other sublevel structs that has their own state and auxiliary variables are stored as fields
# To avoid any confusion, it is recommended to use only state and auxiliary variables in the main struct, or use only the fields of sublevel structs

module Namespace

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using NetcdfIO: read_nc

using ..EmeraldIO.Text: read_csv
using ..EmeraldPhysics.Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, GRAVITY, M_H₂O, P_ATM, T₀, T₂₅, ρ_H₂O


# Please do not use V1/V2/V3 files here as they do not contain the Phi_PSI and Phi_PSII variables
const LAND_2017     = artifact"land_model_spectrum_V4" * "/clima_land_spectra_2017.nc";
const LAND_2021     = artifact"land_model_spectrum_V4" * "/clima_land_spectra_2021.nc";
const LAND_2017_1NM = artifact"land_model_spectrum_V4" * "/clima_land_spectra_1nm_2017.nc";
const LAND_2021_1NM = artifact"land_model_spectrum_V4" * "/clima_land_spectra_1nm_2021.nc";
const SOIL_TEXT     = read_csv("$(@__DIR__)/../../data/SOIL-TEXTURE.csv");


# The configuration of the SPAC system
include("config/spectra.jl");
include("config/trace.jl");

include("config/config.jl");


# General methods (for users to choose from)
include("method/beta.jl");
include("method/colimit.jl");
include("method/kinetics.jl");
include("method/mode.jl");
include("method/pv.jl");
include("method/soil.jl");
include("method/stomata.jl");
include("method/xylem.jl");


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
include("leaf/fluorescence.jl");
include("leaf/flux.jl");

include("leaf/photosynthesis.jl");

include("leaf/leaf.jl");


# Environment
include("environment/air.jl");
include("environment/radiation.jl");

include("environment/meteorology.jl");






include("geometry.jl");

include("canopy.jl");

include("spac.jl");


end
