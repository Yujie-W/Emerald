module SPAC

using Statistics: mean

using ..EmeraldPhysics.Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, GAS_R, M_H₂O, P_ATM, T₀, Λ_THERMAL_H₂O, ρ_H₂O, ρg_MPa
using ..EmeraldPhysics.Optics: photon

using ..CanopyOptics: MODIS_EVI, MODIS_NDVI, MODIS_NIRv, OCO2_SIF759, OCO2_SIF770, TROPOMI_SIF683, TROPOMI_SIF740
using ..CanopyOptics: canopy_fluorescence!, canopy_radiation!, longwave_radiation!, soil_albedo!
using ..EnergyBudget: heat_capacitance, spac_energy_budget!, spac_energy_flow!
using ..LeafOptics: plant_leaf_spectra!
using ..Namespace: XylemHydraulicsAuxilNSS
using ..Namespace: JunctionCapacitor, Leaf, Root, SoilBulk, SoilLayer, Stem
using ..Namespace: AirLayer, GCO₂Mode, MultiLayerSPAC, SPACConfiguration, MultiLayerSPACState
using ..Photosynthesis: plant_photosynthesis!
using ..PhysicalChemistry: latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using ..PlantHydraulics: capacitance_pressure, flow_out, plant_flow_profile!, plant_pressure_profile!, plant_water_budget!, read_β
using ..SoilHydraulics: relative_soil_k, soil_budgets!, soil_profiles!, soil_ψ_25
using ..StomatalModels: stomatal_conductance!, stomatal_conductance_profile!


include("instructions/initialize.jl");
include("instructions/update_auxil.jl");

include("budget.jl");
include("initialize.jl");
include("model.jl");
include("quantity.jl");
include("state.jl");
include("update.jl");


end; # module
