module SPAC

using Statistics: mean

using ..Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, P_ATM, T₀, ρ_H₂O
using ..PhysicalChemistry: latent_heat_vapor, saturation_vapor_pressure
using ..Namespace: AirLayer, GCO₂Mode, MultiLayerSPAC, SPACConfiguration, MultiLayerSPACState
using ..LeafOptics: leaf_spectra!
using ..CanopyOptics: MODIS_EVI, MODIS_NDVI, MODIS_NIRv, OCO2_SIF759, OCO2_SIF770, TROPOMI_SIF683, TROPOMI_SIF740
using ..CanopyOptics: canopy_fluorescence!, canopy_radiation!, longwave_radiation!, soil_albedo!
using ..Photosynthesis: leaf_photosynthesis!
using ..SoilHydraulics: soil_budget!
using ..PlantHydraulics: flow_out, plant_energy!, xylem_flow_profile!, xylem_pressure_profile!, β_factor
using ..StomatalModels: stomatal_conductance!, stomatal_conductance_profile!


include("spac/budget.jl");
include("spac/initialize.jl");
include("spac/model.jl");
include("spac/quantity.jl");
include("spac/state.jl");
include("spac/update.jl");


end # module
