module SPAC

using Statistics: mean

using ..EmeraldPhysics.Constant: CP_D_MOL, CP_L, CP_L_MOL, CP_V_MOL, F_N₂, F_O₂, GAS_R, M_H₂O, P_ATM, T₀, Λ_THERMAL_H₂O, ρ_H₂O, ρg_MPa
using ..EmeraldPhysics.Optics: photon

using ..CanopyOptics: canopy_radiation!, fluorescence_spectrum!, longwave_radiation!, reflection_spectrum!, sensor_geometry!, soil_albedo!
using ..EnergyBudget: heat_capacitance, spac_energy_budget!, spac_energy_flow!
using ..LeafOptics: plant_leaf_spectra!
using ..Namespace: XylemHydraulicsAuxilNSS
using ..Namespace: JunctionCapacitor, Leaf, Root, SoilBulk, SoilLayer, Stem, ReferenceSpectra
using ..Namespace: AirLayer, GCO₂Mode, MultiLayerCanopy, BulkSPAC, SPACConfiguration, MultiLayerSPACState
using ..Photosynthesis: plant_photosynthesis!
using ..PhysicalChemistry: latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using ..PlantHydraulics: capacitance_pressure, flow_out, plant_flow_profile!, plant_pressure_profile!, plant_water_budget!
using ..SoilHydraulics: relative_soil_k, soil_budgets!, soil_profiles!, soil_ψ_25
using ..StomatalModels: limit_stomatal_conductance!, stomatal_conductance!, stomatal_conductance_profile!, read_β


# general instructions to run SPAC
include("instructions/initialize.jl");
include("instructions/prescribe.jl");
include("instructions/update_auxil.jl");


# time stepper
include("timestepper/timer.jl");

include("timestepper/stepper.jl");


# quantities of SPAC
include("quantity/beta.jl");
include("quantity/et.jl");
include("quantity/etr.jl");
include("quantity/gpp.jl");
include("quantity/npp.jl");
include("quantity/par.jl");
include("quantity/sif.jl");
include("quantity/yield.jl");

include("quantity/spectrum.jl");

include("quantity/modis.jl");
include("quantity/oco.jl");
include("quantity/tropomi.jl");


include("model.jl");
include("state.jl");


end; # module
