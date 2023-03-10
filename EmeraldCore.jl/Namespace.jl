module NameSpace

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS

using ..Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, M_H₂O, P_ATM, T₀, T₂₅, ρ_H₂O, ρg_MPa
using ..EmeraldIO.Netcdf: read_nc

const LAND_2017 = artifact"land_model_spectrum_V2" * "/clima_land_spectra_2017.nc";
const LAND_2021 = artifact"land_model_spectrum_V2" * "/clima_land_spectra_2021.nc";


include("namespace/air.jl"        )
include("namespace/colimit.jl"    )
include("namespace/meteorology.jl")
include("namespace/soil.jl"       )
include("namespace/trace.jl"      )


end
