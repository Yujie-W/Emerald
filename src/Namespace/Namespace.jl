module Namespace

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using LazyArtifacts
using PkgUtility.DataIO: read_csv, read_jld2, save_jld2!
using PkgUtility.MathTools: NewtonBisectionMethod, SolutionTolerance, find_zero, interpolate_data
using PkgUtility.RecursiveTools: sync_struct!
using PkgUtility.UniversalConstants: TraceGasAir, TraceGasCH₄, TraceGasCO₂, TraceGasH₂O, TraceGasN₂, TraceGasO₂, TraceLiquidH₂O
using PkgUtility.UniversalConstants: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, GRAVITY, M_H₂O, P_ATM, T₀, T₂₅, ρ_H₂O


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
const SOIL_TEXTURE     = read_csv("$(@__DIR__)/../../data/SOIL-TEXTURE.csv");


# Configurations
include("config/config-info.jl");
include("config/constants/photosynthesis-rate-constant.jl");
include("config/constants/reference-spectra.jl");
include("config/constants.jl");
include("config/dimensions.jl");
include("config/features.jl");
include("config/methods/colimitation-method.jl");
include("config/methods/colimitation-method-settings.jl");
include("config/methods/fluorescence-model.jl");
include("config/methods/fluorescence-model-settings.jl");
include("config/methods/fluorescence-spectra.jl");
include("config/methods/photosynthesis-model.jl");
include("config/methods/soil-albedo.jl");
include("config/methods/stomatal-model-beta.jl");
include("config/methods/stomatal-model.jl");
include("config/methods/temperature-dependency.jl");
include("config/methods/temperature-dependency-settings.jl");
include("config/methods.jl");
include("config.jl");

# General methods (for users to choose from)
include("method/lidf.jl");
include("method/pv.jl");
include("method/soil.jl");
include("method/xylem.jl");

# Plant hydraulics (dependent on config and method)
include("plant/xylem/energy.jl");
include("plant/xylem/junction.jl");
include("plant/xylem/xylem.jl");

# Soil
include("soil/bulk.jl");
include("soil/layer.jl");

# Root system (dependent on xylem)
include("plant/root/rhizosphere.jl");
include("plant/root/root.jl");

# Stem system (dependent on xylem)
include("plant/stem/stem.jl");

# Leaf system (dependent on xylem)
include("plant/leaf/biophysics.jl");
include("plant/leaf/energy.jl");
include("plant/leaf/extraxylem.jl");
include("plant/leaf/leafflux.jl");
include("plant/leaf/layerflux.jl");
include("plant/leaf/photosynthesis.jl");
include("plant/leaf/layer.jl");
include("plant/leaf/leaf.jl");

# Canopy
include("plant/canopy/clumping.jl");
include("plant/canopy/sensor_geometry.jl");
include("plant/canopy/structure.jl");
include("plant/canopy/sun_geometry.jl");
include("plant/canopy/canopy.jl");

# Environment
include("environment/air.jl");
include("environment/radiation.jl");
include("environment/meteorology.jl");

# SPAC
include("spac/cache.jl");
include("spac/info.jl");
include("spac/memory.jl");
include("spac/pool.jl");
include("spac/plant.jl");
include("spac/bulk.jl");
include("general.jl");


end # module Namespace
