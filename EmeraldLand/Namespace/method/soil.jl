# This file contains soil types

#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2021-Sep-30: add abstract type for soil vulnerability curve
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractSoilVC:
- [`BrooksCorey`](@ref)
- [`VanGenuchten`](@ref)

"""
abstract type AbstractSoilVC{FT<:AbstractFloat} end;


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Sep-30: define this structure with no default constructor
#     2022-Oct-19: make struct mutable
#     2023-May-12: add BrooksCorey constructor
# Sources
#     https://ral.ucar.edu/sites/default/files/public/product-tool/unified-noah-lsm/parameters/SOILPARM.TBL
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Brooks Corey soil parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct BrooksCorey{FT<:AbstractFloat} <:AbstractSoilVC{FT}
    # General model information
    "Soil b"
    B::FT
    "Maximum soil hydraulic conductivity at 25 °C `[mol m⁻¹ s⁻¹ MPa⁻¹]`"
    K_MAX::FT
    "Soil type"
    TYPE::String
    "Potential at saturation `[MPa]`"
    Ψ_SAT::FT
    "Saturated soil volumetric water content"
    Θ_SAT::FT
    "Residual soil volumetric water content"
    Θ_RES::FT
end;

"""

    BrooksCorey{FT}(catg::Int) where {FT}

Return a BrooksCorey soil VC, given
- `catg` Soil texture catergory (must be within [1,19])

"""
BrooksCorey{FT}(catg::Int) where {FT} = (
    @assert 1 <= catg <= 19 "Soil texture catergory must be within 1 to 19!";

    return BrooksCorey{FT}(
                B     = SOIL_TEXT.BB[catg],
                K_MAX = SOIL_TEXT.SATDK[catg] / GRAVITY(FT) * 1e6 / M_H₂O(FT),
                TYPE  = SOIL_TEXT.NAME[catg],
                Ψ_SAT = SOIL_TEXT.SATPSI[catg] * ρ_H₂O(FT) * GRAVITY(FT) * 1e-6,
                Θ_SAT = SOIL_TEXT.MAXSMC[catg],
                Θ_RES = SOIL_TEXT.REFSMC[catg]
    )
);

sync_state!(state_from::BrooksCorey{FT}, state_to::BrooksCorey{FT}) where {FT} = (
    state_to.B = state_from.B;
    state_to.K_MAX = state_from.K_MAX;
    state_to.TYPE = state_from.TYPE;
    state_to.Ψ_SAT = state_from.Ψ_SAT;
    state_to.Θ_SAT = state_from.Θ_SAT;
    state_to.Θ_RES = state_from.Θ_RES;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2021-Sep-30: define this structure with two default constructors from an incomplete parameter set
#     2021-Sep-30: add constructor function
#     2022-Jul-13: remove a constructor method
#     2022-Oct-19: make struct mutable
#     2023-Apr-08: add another dataset for van Genuchten parameters
# Sources
#     Cosby et al. (1984) A statistical exploration of the relationships of soil moisture characteristics to the physical properties of soils
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

van Genuchten soil parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct VanGenuchten{FT<:AbstractFloat} <:AbstractSoilVC{FT}
    # General model information
    "Maximum soil hydraulic conductivity at 25 °C `[mol m⁻¹ s⁻¹ MPa⁻¹]`"
    K_MAX::FT
    "Soil n is Measure of the pore-size distribution"
    N::FT
    "Soil type"
    TYPE::String
    "Soil α is related to the inverse of the air entry suction, α > 0"
    α::FT
    "Residual soil volumetric water content"
    Θ_RES::FT
    "Saturated soil volumetric water content"
    Θ_SAT::FT

    # Parameters based on the ones above
    "Soil m = 1 - 1/n"
    M::FT = 1 - 1 / N
end;


"""

    VanGenuchten{FT}(name::String) where {FT}

Constructor for [`VanGenuchten`](@ref), given
- `name` Soil type name, need to be Sand, Loamy Sand, Sandy Loam, Loam (default), Sandy Clay Loam, Silt Loam, Silt, Clay Loam, Silty Clay Loam, Sandy Clay, Silty Clay, and Clay.

"""
VanGenuchten{FT}(name::String) where {FT} = (
    #=
    # Parameters from Loam soil
    _p = [ 367.3476, 1.56, 0.43, 0.078, exp(-0.32) * 0.0254 / 3600];

    # switch name
    if name=="Sand"
        _p = [1479.5945, 2.68, 0.43, 0.045, exp( 0.82) * 0.0254 / 3600];
    elseif name=="Loamy Sand"
        _p = [1265.3084, 2.28, 0.41, 0.057, exp( 0.30) * 0.0254 / 3600];
    elseif name=="Sandy Loam"
        _p = [ 765.3075, 1.89, 0.41, 0.065, exp(-0.13) * 0.0254 / 3600];
    elseif name=="Loam"
        _p = [ 367.3476, 1.56, 0.43, 0.078, exp(-0.32) * 0.0254 / 3600];
    elseif name=="Sandy Clay Loam"
        _p = [ 602.0419, 1.48, 0.39, 0.100, exp(-0.20) * 0.0254 / 3600];
    elseif name=="Silt Loam"
        _p = [ 204.0820, 1.41, 0.45, 0.067, exp(-0.40) * 0.0254 / 3600];
    elseif name=="Silt"
        _p = [ 163.2656, 1.37, 0.46, 0.034, exp(-0.63) * 0.0254 / 3600];   # this k is guessed, must be within [-0.72, -0.54]
    elseif name=="Clay Loam"
        _p = [ 193.8779, 1.31, 0.41, 0.095, exp(-0.40) * 0.0254 / 3600];
    elseif name=="Silty Clay Loam"
        _p = [ 102.0410, 1.23, 0.43, 0.089, exp(-0.54) * 0.0254 / 3600];
    elseif name== "Sandy Clay"
        _p = [ 275.5107, 1.23, 0.38, 0.100, exp( 0.01) * 0.0254 / 3600];
    elseif name=="Silty Clay"
        _p = [  51.0205, 1.09, 0.36, 0.070, exp(-0.72) * 0.0254 / 3600];
    elseif name=="Clay"
        _p = [  81.6328, 1.09, 0.38, 0.068, exp(-0.86) * 0.0254 / 3600];    # K from Light clay in Cosby et al. (1984)
    else
        @warn "Soil type $(name) not recognized, use Loam instead.";
        name = "Loam";
    end;
    =#
    # https://structx.com/Soil_Properties_007.html
    # Parameters from Loam soil
    p = [ 367.3476, 1.56, 0.43, 0.078, 7.19e-6];

    # switch name
    if name=="Sand"
        p = [1479.5945, 2.68, 0.43, 0.045, 1.76e-4];
    elseif name=="Loamy Sand"
        p = [1265.3084, 2.28, 0.41, 0.057, 1.56e-4];
    elseif name=="Sandy Loam"
        p = [ 765.3075, 1.89, 0.41, 0.065, 3.45e-5];
    elseif name=="Loam"
        p = [ 367.3476, 1.56, 0.43, 0.078, 6.94e-6];
    elseif name=="Sandy Clay Loam"
        p = [ 602.0419, 1.48, 0.39, 0.100, 6.31e-6];
    elseif name=="Silt Loam"
        p = [ 204.0820, 1.41, 0.45, 0.067, 7.19e-6];
    elseif name=="Silt"
        p = [ 163.2656, 1.37, 0.46, 0.034, 7.19e-6];
    elseif name=="Clay Loam"
        p = [ 193.8779, 1.31, 0.41, 0.095, 2.45e-6];
    elseif name=="Silty Clay Loam"
        p = [ 102.0410, 1.23, 0.43, 0.089, 1.70e-6];
    elseif name== "Sandy Clay"
        p = [ 275.5107, 1.23, 0.38, 0.100, 2.17e-6];
    elseif name=="Silty Clay"
        p = [  51.0205, 1.09, 0.36, 0.070, 1.02e-6];
    elseif name=="Clay"
        p = [  81.6328, 1.09, 0.38, 0.068, 1.28e-6];
    else
        @warn "Soil type $(name) not recognized, use Loam instead.";
        name = "Loam";
    end;

    # return a new struct
    return VanGenuchten{FT}(K_MAX = p[5] / GRAVITY(FT) * 1e6 / M_H₂O(FT), N = p[2], TYPE = name, α = p[1], Θ_RES = p[4], Θ_SAT = p[3])
);

sync_state!(state_from::VanGenuchten{FT}, state_to::VanGenuchten{FT}) where {FT} = (
    state_to.K_MAX = state_from.K_MAX;
    state_to.N = state_from.N;
    state_to.TYPE = state_from.TYPE;
    state_to.α = state_from.α;
    state_to.Θ_RES = state_from.Θ_RES;
    state_to.Θ_SAT = state_from.Θ_SAT;
    state_to.M = state_from.M;

    return nothing
);
