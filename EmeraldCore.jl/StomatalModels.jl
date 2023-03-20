module StomatalModels

using ..Constant: CP_D_MOL, K_STEFAN, M_H₂O
using ..Namespace: AbstractSoilVC, AbstractXylemVC
using ..Namespace: AbstractStomataModel, AndereggSM, BallBerrySM, BetaParameterG1, BetaParameterVcmax, EllerSM, GentineSM, LeuningSM, MedlynSM, SperrySM, WangSM, Wang2SM
using ..Namespace: AirLayer, BroadbandLeafBiophysics, C4VJPModel, GCO₂Mode, HyperspectralLeafBiophysics, Leaf, LeafHydraulics, Leaves1D, Leaves2D
using ..Namespace: MonoElementSPAC, MonoMLTreeSPAC
using ..Photosynthesis: leaf_photosynthesis!, ∂R∂T
using ..PhysicalChemistry: latent_heat_vapor, relative_diffusive_coefficient, relative_surface_tension, relative_viscosity, saturation_vapor_pressure
using ..PlantHydraulics: relative_hydraulic_conductance, ∂E∂P
using ..SoilHydraulics: relative_hydraulic_conductance


include("stomata/conductance.jl")
include("stomata/empirical.jl"  )
include("stomata/limits.jl"     )
include("stomata/optimality.jl" )


end # module
