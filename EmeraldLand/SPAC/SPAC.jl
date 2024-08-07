module SPAC

using Statistics: mean

using ..EmeraldPhysics.Constant: F_N₂, F_O₂, GAS_R, T₀, Λ_THERMAL_H₂O
using ..EmeraldPhysics.Optics: photon

using ..CanopyOptics: canopy_structure!, canopy_structure_aux!, longwave_radiation!, shortwave_radiation!, soil_albedo!, sun_geometry!, sun_geometry_aux!
using ..CanopyOptics: fluorescence_spectrum!, reflection_spectrum!, sensor_geometry!, sensor_geometry_aux!
using ..EnergyBudget: heat_capacitance, spac_energy_budget!, spac_energy_flow!
using ..LeafOptics: plant_leaf_spectra!
using ..Namespace: ReferenceSpectra, ShortwaveRadiation
using ..Namespace: AcMethodC3VcmaxPi, AcMethodC4Vcmax
using ..Namespace: AjMethodC3JmaxPi, AjMethodC3VqmaxPi, AjMethodC4JPSII
using ..Namespace: ApMethodC3Inf, ApMethodC3Vcmax, ApMethodC4VcmaxPi, ApMethodC4VpmaxPi
using ..Namespace: GeneralC3Trait, GeneralC4Trait
using ..Namespace: ExtraXylemCapacitorState, XylemHydraulicsAuxilNSS, XylemHydraulicsTrait, LeafBioTrait, LeafEnergyState, LeafEnergySDAuxil
using ..Namespace: CanopyLayer, JunctionCapacitor, Leaf, Plant, Root, Stem
using ..Namespace: CanopyStructure, CanopyStructureTrait, CanopyStructureTDAuxil, MultiLayerCanopy
using ..Namespace: AirLayer, AirLayerState, AirLayerSDAuxil, AirLayerTDAuxil
using ..Namespace: SoilBulk, SoilLayer, SoilLayerState, SoilLayerSDAuxil, SoilLayerTrait, SoilLayerTDAuxil
using ..Namespace: BulkSPAC, BulkSPACStates, SPACCache, SPACConfiguration, sync_state!
using ..Photosynthesis: plant_photosynthesis!
using ..PhysicalChemistry: relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using ..PlantHydraulics: capacitance_pressure, flow_out, plant_flow_profile!, plant_pressure_profile!, plant_water_budget!, xylem_conductance, xylem_pressure
using ..SoilHydraulics: relative_soil_k, soil_budgets!, soil_profiles!, soil_ψ_25
using ..StomatalModels: limit_stomatal_conductance!, read_β, stomatal_conductance!, stomatal_conductance_profile!, β_factor!


# general instructions to run SPAC
include("instructions/aux_dull.jl");
include("instructions/aux_state.jl");
include("instructions/aux_step.jl");
include("instructions/aux_substep.jl");
include("instructions/aux_trait.jl");
include("instructions/initialize.jl");
include("instructions/model_step.jl");
include("instructions/model_substep.jl");
include("instructions/prescribe_air.jl");
include("instructions/prescribe_plant.jl");
include("instructions/prescribe_soil.jl");
include("instructions/push_t_history.jl");
include("instructions/shed_leaves.jl");


# time stepper
include("timestepper/model.jl");
include("timestepper/stepper.jl");
include("timestepper/timer.jl");


# quantities of SPAC
include("quantity/beta.jl");
include("quantity/et.jl");
include("quantity/gpp.jl");
include("quantity/goes.jl");
include("quantity/hydraulics.jl");
include("quantity/npp.jl");
include("quantity/par.jl");
include("quantity/sif.jl");
include("quantity/yield.jl");

include("quantity/modis.jl");
include("quantity/oco.jl");
include("quantity/spectrum.jl");
include("quantity/tropomi.jl");


end; # module
