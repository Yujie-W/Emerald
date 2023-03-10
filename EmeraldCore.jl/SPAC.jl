module SPAC

using Statistics: mean

using ..Constant: CP_L, CP_L_MOL, T₀, ρ_H₂O
using ..PhysicalChemistry: saturation_vapor_pressure
using ..Namespace: AirLayer, GCO₂Mode, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
using ..LeafOptics: leaf_spectra!
using ..CanopyOptics: canopy_fluorescence!, canopy_radiation!, soil_albedo!
using ..Photosynthesis: leaf_photosynthesis!
using ..Soil: soil_budget!
using ..PlantHydraulics: flow_out, plant_energy!, xylem_flow_profile!, xylem_pressure_profile!, β_factor
using ..StomatalModels: stomatal_conductance!


include("spac/budget.jl"    )
include("spac/initialize.jl")
include("spac/model.jl"     )
include("spac/quantity.jl"  )
include("spac/update.jl"    )


end # module
