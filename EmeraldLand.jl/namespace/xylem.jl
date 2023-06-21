#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Feb-01: add abstract type for vulnerability curve
#     2022-Apr-20: rename type to AbstractXylemVC
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractXylemVC:
- [`LogisticVC`](@ref)
- [`PowerVC`](@ref)
- [`WeibullVC`](@ref)
- [`ComplexVC`](@ref)

"""
abstract type AbstractXylemVC{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Feb-01: add LogisticVC (use mutable so as to fit)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Modified logistic function for vulnerability curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LogisticVC{FT<:AbstractFloat} <: AbstractXylemVC{FT}
    # General model information
    "Multiplier to exponential component"
    A::FT = 1
    "Multiplier to pressure `[MPa⁻¹]`"
    B::FT = 1
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Feb-01: add PowerVC (use mutable so as to fit)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Power function for vulnerability curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct PowerVC{FT<:AbstractFloat} <: AbstractXylemVC{FT}
    # General model information
    "Multiplier to power component `[MPa⁻ᵇ]`"
    A::FT = 1
    "Power to pressure"
    B::FT = 1
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Feb-01: add WeibullVC (use mutable so as to fit)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Weibull cumulative distribution function for vulnerability curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct WeibullVC{FT<:AbstractFloat} <: AbstractXylemVC{FT}
    # General model information
    "Numerator in the exponential component `[MPa]`"
    B::FT = 2
    "Power to pressure component"
    C::FT = 5
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Feb-01: add ComplexVC (use mutable so as to fit)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

A complex struct for segmented vulnerability curve such as dual Weibull function

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ComplexVC{FT<:AbstractFloat} <: AbstractXylemVC{FT}
    # General model information
    "Percentages of each VC component"
    PS::Vector{FT} = FT[0.5, 0.5]
    "Vector of vulnerability curve components"
    VCS::Union{Vector{LogisticVC{FT}}, Vector{PowerVC{FT}}, Vector{WeibullVC{FT}}, Vector{AbstractXylemVC{FT}}} = WeibullVC{FT}[WeibullVC{FT}(), WeibullVC{FT}(B = 3)]
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Apr-20: add abstract type for pressure volume curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractPVCurve:
- [`LinearPVCurve`](@ref)
- [`SegmentedPVCurve`](@ref)

"""
abstract type AbstractPVCurve{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Apr-20: add linear PV curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains information for linear PV curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LinearPVCurve{FT<:AbstractFloat} <: AbstractPVCurve{FT}
    # General model information
    "Conductance for refilling (relative to maximum) `[MPa⁻¹ s⁻¹]`"
    K_REFILL::FT = 1e4
    "Slope of the linear PV curve (relative to maximum) `[MPa⁻¹]`"
    SLOPE::FT = 0.2
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-24: add segmented PV curve
#     2022-Jul-20: rename Ε_BULK (greek) to ϵ_BULK
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains information for segmented PV curve

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SegmentedPVCurve{FT<:AbstractFloat} <: AbstractPVCurve{FT}
    # General model information
    "n_o / maximum V `[mol m⁻³]`"
    C_ALL::FT = 300
    "Conductance for refilling (relative to maximum) `[MPa⁻¹ s⁻¹]`"
    K_REFILL::FT = 1e-4
    "Apoplastic water content relative to maximum water volume"
    RWC_APO::FT = 0.2
    "Relative water content at turgor loss point"
    RWC_TLP::FT = 0.8
    "Bulk modulus of elasticity `[MPa]`"
    ϵ_BULK::FT = 20
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-27: add abstract type for steady and non-steady flow
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractFlowProfile:
- [`NonSteadyStateFlow`](@ref)
- [`SteadyStateFlow`](@ref)

"""
abstract type AbstractFlowProfile{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-27: add non-steady state flow struct
#     2022-May-31: remove unneeded variables from struct
#     2022-May-31: add field: f_sum
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains stem hydraulic system flow rates at non-steady state

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct NonSteadyStateFlow{FT<:AbstractFloat} <: AbstractFlowProfile{FT}
    # Dimensions
    "Dimension of capaciatance elements"
    DIM_CAPACITY::Int = 1

    # Diagnostic variables
    "Flow rate into the organ `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    f_in::FT = 0
    "Flow rate out of the organ `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    f_out::FT = 0

    # Cache variables
    "Vector of xylem water flow `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    _f_element::Vector{FT} = zeros(FT, DIM_CAPACITY)
    "Vector of buffer water flow `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    _f_buffer::Vector{FT} = zeros(FT, DIM_CAPACITY)
    "Vector of sum buffer water flow `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    _f_sum::Vector{FT} = zeros(FT, DIM_CAPACITY)
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-27: add non-steady state flow struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains stem hydraulic system flow rates at steady state

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SteadyStateFlow{FT<:AbstractFloat} <: AbstractFlowProfile{FT}
    # Diagnostic variables
    "Flow rate through the organ `[mol s⁻¹]` (for root and stem) or `[mol m⁻² s⁻¹]` (for leaf)"
    flow::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-May-24: add abstract type for hydraulic organ
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractHydraulicSystem:
- [`LeafHydraulics`](@ref)
- [`RootHydraulics`](@ref)
- [`StemHydraulics`](@ref)

"""
abstract type AbstractHydraulicSystem{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: add leaf hydraulic system
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-07: add e_crit as a field
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains leaf hydraulic system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafHydraulics{FT<:AbstractFloat} <: AbstractHydraulicSystem{FT}
    # Dimensions
    "Dimension of xylem slices"
    DIM_XYLEM::Int = 5

    # General information of the hydraulic system
    "Leaf area `[m²]`"
    AREA::FT = 1500
    "Maximal extra-xylary hydraulic conductance per leaf area `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_OX::FT = 100
    "Maximal leaf xylem hydraulic conductance per leaf area `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_SLA::FT = 0.04
    "Total capaciatance at Ψ = 0 `[mol m⁻²]`"
    V_MAXIMUM::FT = 20

    # Embedded structures
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}} = SteadyStateFlow{FT}()
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}} = SegmentedPVCurve{FT}()
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}} = WeibullVC{FT}()

    # Prognostic variables (used for ∂y∂t)
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Current capaciatance at Ψ_leaf `[mol m⁻²]`"
    v_storage::FT = V_MAXIMUM

    # Prognostic variables (not used for ∂y∂t)
    "Leaf end water pressure at air-water interface `[MPa]`"
    p_leaf::FT = 0
    "Leaf xylem water pressure at the leaf base (upstream) `[MPa]`"
    p_ups::FT = 0

    # Cache variables
    "Critical flow rate `[mol s⁻¹ m⁻²]`"
    _e_crit::FT = 0
    "Vector of leaf kr history per element `[-]`"
    _k_history::Vector{FT} = ones(FT, DIM_XYLEM)
    "Leaf xylem water pressure at the downstream end of leaf xylem `[MPa]`"
    _p_dos::FT = 0
    "Vector of xylem water pressure `[MPa]`"
    _p_element::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Pressure of storage"
    _p_storage::FT = 0
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: add root hydraulic system
#     2022-May-25: rename PV to PVC to be consistent with LeafHydraulics
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-19: add dimension control to struct
#     2022-Oct-20: remove field SH (use that in SOIL in the future)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains root hydraulic system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct RootHydraulics{FT<:AbstractFloat} <: AbstractHydraulicSystem{FT}
    # Dimensions
    "Dimension of xylem slices"
    DIM_XYLEM::Int = 5

    # General information of the hydraulic system
    "Root cross-section area `[m²]`"
    AREA::FT = 1
    "Rhizosphere  conductance `[mol s⁻¹ MPa⁻¹]`"
    K_RHIZ::FT = 5e14
    "Maximal xylem hydraulic conductivity `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_X::FT = 25
    "Length `[m]`"
    L::FT = 1
    "Maximal storage per element `[mol]`"
    V_MAXIMUM::Vector{FT} = AREA * L / DIM_XYLEM * 6000 * ones(FT, DIM_XYLEM)
    "Root z difference `[m]`"
    ΔH::FT = 1

    # Embedded structures
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}} = SteadyStateFlow{FT}()
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}} = LinearPVCurve{FT}()
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}} = WeibullVC{FT}()

    # Prognostic variables (used for ∂y∂t)
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Storage per element `[mol]`"
    v_storage::Vector{FT} = V_MAXIMUM .* 1

    # Prognostic variables (not used for ∂y∂t)
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos::FT = 0
    "Soil matrix potential `[MPa]`"
    p_ups::FT = 0
    "Soil osmotic potential at 298.15 K `[MPa]"
    ψ_osm::FT = 0

    # Cache variables
    "Vector of leaf kr history per element"
    _k_history::Vector{FT} = ones(FT, DIM_XYLEM)
    "Vector of xylem water pressure `[MPa]`"
    _p_element::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Xylem-rhizosphere interface water pressure `[MPa]`"
    _p_rhiz::FT = 0
    "Pressure of storage per element `[MPa]`"
    _p_storage::Vector{FT} = zeros(FT, DIM_XYLEM)
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-May-25: add stem hydraulic system
#     2022-May-25: rename PV to PVC to be consistent with LeafHydraulics
#     2022-May-27: move flow rates to a field FLOW
#     2022-Jun-13: use Union instead of Abstract... for type definition
#     2022-Jul-19: add dimension control to struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct that contains stem hydraulic system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct StemHydraulics{FT<:AbstractFloat} <: AbstractHydraulicSystem{FT}
    # Dimensions
    "Dimension of xylem slices"
    DIM_XYLEM::Int = 5

    # General information of the hydraulic system
    "Xylem cross-section area `[m²]`"
    AREA::FT = 1
    "Maximal xylem hydraulic conductivity (per xylem area per length) `[mol s⁻¹ MPa⁻¹ m⁻²]`"
    K_X::FT = 25
    "Length `[m]`"
    L::FT = 1
    "Maximal storage per element `[mol]`"
    V_MAXIMUM::Vector{FT} = AREA * L / DIM_XYLEM * 6000 * ones(FT, DIM_XYLEM)
    "Root z difference `[m]`"
    ΔH::FT = 1

    # Embedded structures
    "Flow profile"
    FLOW::Union{SteadyStateFlow{FT}, NonSteadyStateFlow{FT}} = SteadyStateFlow{FT}()
    "Pressure volume curve for storage"
    PVC::Union{LinearPVCurve{FT}, SegmentedPVCurve{FT}} = LinearPVCurve{FT}()
    "Vulnerability curve"
    VC::Union{LogisticVC{FT}, PowerVC{FT}, WeibullVC{FT}, ComplexVC{FT}} = WeibullVC{FT}()

    # Prognostic variables (used for ∂y∂t)
    "Vector of xylem water pressure history (normalized to 298.15 K) `[MPa]`"
    p_history::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Storage per element `[mol]`"
    v_storage::Vector{FT} = V_MAXIMUM .* 1

    # Prognostic variables (not used for ∂y∂t)
    "Xylem water pressure at the downstream end of xylem `[MPa]`"
    p_dos::FT = 0
    "Soil matrix potential `[MPa]`"
    p_ups::FT = 0

    # Cache variables
    "Vector of leaf kr history per element"
    _k_history::Vector{FT} = ones(FT, DIM_XYLEM)
    "Vector of xylem water pressure `[MPa]`"
    _p_element::Vector{FT} = zeros(FT, DIM_XYLEM)
    "Pressure of storage per element"
    _p_storage::Vector{FT} = zeros(FT, DIM_XYLEM)
end
