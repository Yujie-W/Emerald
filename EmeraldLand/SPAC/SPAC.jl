module SPAC

using Statistics: mean

using ..EmeraldPhysics.Constant: F_N₂, F_O₂, GAS_R, T₀, Λ_THERMAL_H₂O
using ..EmeraldPhysics.Optics: photon

using ..CanopyOptics: canopy_structure!, inclination_angles!, longwave_radiation!, shortwave_radiation!, soil_albedo!, sun_geometry!
using ..CanopyOptics: fluorescence_spectrum!, reflection_spectrum!, sensor_geometry!
using ..EnergyBudget: heat_capacitance, spac_energy_budget!, spac_energy_flow!
using ..LeafOptics: plant_leaf_spectra!
using ..Namespace: ReferenceSpectra
using ..Namespace: GCO₂Mode
using ..Namespace: XylemHydraulicsAuxilNSS
using ..Namespace: JunctionCapacitor, Leaf, Root, Stem
using ..Namespace: MultiLayerCanopy
using ..Namespace: AirLayer, SoilBulk, SoilLayer
using ..Namespace: BulkSPAC, SPACConfiguration
using ..Photosynthesis: plant_photosynthesis!
using ..PhysicalChemistry: relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using ..PlantHydraulics: capacitance_pressure, flow_out, plant_flow_profile!, plant_pressure_profile!, plant_water_budget!
using ..SoilHydraulics: relative_soil_k, soil_budgets!, soil_profiles!, soil_ψ_25
using ..StomatalModels: limit_stomatal_conductance!, read_β, stomatal_conductance!, stomatal_conductance_profile!


# general instructions to run SPAC
include("instructions/initialize.jl");
include("instructions/prescribe.jl");
include("instructions/update_auxil.jl");


# time stepper
include("timestepper/timer.jl");

include("timestepper/stepper.jl");

include("timestepper/model.jl");


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


# This file is to be refactored when reaching the EmeraldEarth stage
using.. Namespace: MultiLayerSPACState
include("state.jl");


end; # module
