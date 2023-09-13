module Namespace

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using NetcdfIO: read_nc

using ..EmeraldIO.Text: read_csv

using ..Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, GRAVITY, M_H₂O, P_ATM, T₀, T₂₅, ρ_H₂O


# Please do not use V1/V2/V3 files here as they do not contain the Phi_PSI and Phi_PSII variables
const LAND_2017     = artifact"land_model_spectrum_V4" * "/clima_land_spectra_2017.nc";
const LAND_2021     = artifact"land_model_spectrum_V4" * "/clima_land_spectra_2021.nc";
const LAND_2017_1NM = artifact"land_model_spectrum_V4" * "/clima_land_spectra_1nm_2017.nc";
const LAND_2021_1NM = artifact"land_model_spectrum_V4" * "/clima_land_spectra_1nm_2021.nc";
const SOIL_TEXT     = read_csv("$(@__DIR__)/../data/SOIL-TEXTURE.csv");


include("namespace/pigment.jl");
include("namespace/radiation.jl");
include("namespace/trace.jl");

include("namespace/config.jl");

include("namespace/air.jl");
include("namespace/colimit.jl");
include("namespace/geometry.jl");
include("namespace/kinetics.jl");
include("namespace/meteorology.jl");
include("namespace/soil.jl");
include("namespace/stomata.jl");
include("namespace/xylem.jl");

include("namespace/canopy.jl");
include("namespace/leaf.jl");
include("namespace/root.jl");
include("namespace/stem.jl");

include("namespace/spac.jl");


end
