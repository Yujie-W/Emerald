#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2024-Jul-24: add cache struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct to save SPAC cache arrays and solvers

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct SPACCache{FT}
    # Cache vectors with the length of Layer
    "Cache vector with the length of Layer"
    cache_layer_1::Vector{FT}

    # Cache vectors with the length of SIFE
    "Cache vector with the length of SIFE"
    cache_sife_1::Vector{FT}
    "Cache vector with the length of SIFE"
    cache_sife_2::Vector{FT}
    "Cache vector with the length of SIFE"
    cache_sife_3::Vector{FT}

    # Cache vectors with the length of wavelength
    "Cache vector with the length of wavelength"
    cache_wl_1::Vector{FT}
    "Cache vector with the length of wavelength"
    cache_wl_2::Vector{FT}
    "Cache vector with the length of wavelength"
    cache_wl_3::Vector{FT}
    "Cache vector with the length of wavelength"
    cache_wl_4::Vector{FT}
    "Cache vector with the length of wavelength"
    cache_wl_5::Vector{FT}

    # Matrix with dims of INCL and AZI
    "Cache matrix with the dims of INCL and AZI"
    cache_incl_azi::Matrix{FT}

    # solvers and tolerances
    "NewtonBisectionMethod solver"
    solver_nb::NewtonBisectionMethod{FT}
    "SolutionTolerance for NewtonBisectionMethod"
    stol_nb::SolutionTolerance{FT}
end;

SPACCache{FT}(dim_azi::Int, dim_incl::Int, dim_layer::Int, dim_sif::Int, dim_sife::Int, dim_wl::Int) where {FT} = SPACCache{FT}(
    cache_layer_1 = zeros(FT, dim_layer),

    cache_sife_1 = zeros(FT, dim_sife),
    cache_sife_2 = zeros(FT, dim_sife),
    cache_sife_3 = zeros(FT, dim_sife),

    cache_wl_1 = zeros(FT, dim_wl),
    cache_wl_2 = zeros(FT, dim_wl),
    cache_wl_3 = zeros(FT, dim_wl),
    cache_wl_4 = zeros(FT, dim_wl),
    cache_wl_5 = zeros(FT, dim_wl),

    cache_incl_azi = zeros(FT, dim_incl, dim_azi),

    solver_nb = NewtonBisectionMethod{FT}(),
    stol_nb   = SolutionTolerance{FT}(eps(FT)*100, 50),
);
