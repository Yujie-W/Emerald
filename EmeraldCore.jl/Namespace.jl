module Namespace

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using StaticArrays: SVector

using ..EmeraldIO.Netcdf: read_nc, size_nc

using ..Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, M_H₂O, P_ATM, T₀, T₂₅, ρ_H₂O, ρg_MPa

const LAND_2017 = artifact"land_model_spectrum_V2" * "/clima_land_spectra_2017.nc";
const LAND_2021 = artifact"land_model_spectrum_V2" * "/clima_land_spectra_2021.nc";


include("namespace/air.jl"      )   # no dims
include("namespace/colimit.jl"  )   # no dims
include("namespace/general.jl"  )   # no dims
include("namespace/geometry.jl" )   # no dims
include("namespace/kinetics.jl" )   # no dims
include("namespace/stomata.jl"  )   # no dims
include("namespace/trace.jl"    )   # no dims

include("namespace/dimension.jl")   # rely on general

include("namespace/canopy.jl"    )  # rely on dims
include("namespace/pigments.jl"  )  # rely on dims
include("namespace/radiation.jl" )  # rely on dims
include("namespace/wavelength.jl")  # rely on dims

include("namespace/config.jl"     ) # rely on pigments, radiation, and wavelength
include("namespace/meteorology.jl") # rely on radiation

include("namespace/soil.jl" )       # rely on config (constructor)
include("namespace/xylem.jl")       # rely on config (constructor)

include("namespace/leaf.jl"  )      # rely on xylem
include("namespace/root.jl"  )      # rely on xylem
include("namespace/stem.jl"  )      # rely on xylem

include("namespace/spac.jl")        # rely on all above


end
