#######################################################################################################################################################################################################
#
# Changes to the struct
# General
#     2024-Jul-24: add cache struct
#     2024-Jul-30: do not bin PPAR if DIM_PPAR_BINS is nothing
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

    # Cache vectors with the length of INCL * AZI + 1
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_1_1::Vector{FT}
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_1_2::Vector{FT}
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_1_3::Vector{FT}
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_1_4::Vector{FT}
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_1_5::Vector{FT}
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_1_6::Vector{FT}
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_2_1::Vector{FT}
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_2_2::Vector{FT}
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_2_3::Vector{FT}
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_2_4::Vector{FT}
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_2_5::Vector{FT}
    "Cache vector with the length of INCL * AZI + 1 or DIM_PPAR_BINS + 1"
    cache_incl_azi_2_6::Vector{FT}


    # Matrix with dims of INCL and AZI
    "Cache matrix with the dims of INCL and AZI"
    cache_incl_azi_1::Matrix{FT}
    "Cache matrix with the dims of INCL and AZI"
    cache_incl_azi_2::Matrix{FT}
    "Cache matrix with the dims of INCL and AZI"
    cache_incl_azi_3::Matrix{FT}
    "Cache matrix with the dims of INCL and AZI"
    cache_incl_azi_4::Matrix{FT}
    "Cache matrix with the dims of INCL and AZI"
    cache_incl_azi_5::Matrix{FT}
    "Cache matrix with the dims of INCL and AZI"
    cache_incl_azi_6::Matrix{FT}
    "Cache matrix with the dims of INCL and AZI"
    cache_incl_azi_7::Matrix{FT}

    # solvers and tolerances
    "NewtonBisectionMethod solver"
    solver_nb::NewtonBisectionMethod{FT}
    "SolutionTolerance for NewtonBisectionMethod"
    stol_nb::SolutionTolerance{FT}
end;

SPACCache{FT}(dim_azi::Int, dim_incl::Int, dim_layer::Int, dim_ppar::Union{Int, Nothing}, dim_sif::Int, dim_sife::Int, dim_wl::Int) where {FT} = (
    cache_dim_ppar = isnothing(dim_ppar) ? dim_incl * dim_azi : dim_ppar;

    return SPACCache{FT}(
                cache_layer_1 = zeros(FT, dim_layer),

                cache_sife_1 = zeros(FT, dim_sife),
                cache_sife_2 = zeros(FT, dim_sife),
                cache_sife_3 = zeros(FT, dim_sife),

                cache_wl_1 = zeros(FT, dim_wl),
                cache_wl_2 = zeros(FT, dim_wl),
                cache_wl_3 = zeros(FT, dim_wl),
                cache_wl_4 = zeros(FT, dim_wl),
                cache_wl_5 = zeros(FT, dim_wl),

                # to used to speed up the computation (PPAR bins)
                cache_incl_azi_1_1 = zeros(FT, cache_dim_ppar+1),
                cache_incl_azi_1_2 = zeros(FT, cache_dim_ppar+1),
                cache_incl_azi_1_3 = zeros(FT, cache_dim_ppar+1),
                cache_incl_azi_1_4 = zeros(FT, cache_dim_ppar+1),
                cache_incl_azi_1_5 = zeros(FT, cache_dim_ppar+1),
                cache_incl_azi_1_6 = zeros(FT, cache_dim_ppar+1),
                cache_incl_azi_2_1 = zeros(FT, cache_dim_ppar+1),
                cache_incl_azi_2_2 = zeros(FT, cache_dim_ppar+1),
                cache_incl_azi_2_3 = zeros(FT, cache_dim_ppar+1),
                cache_incl_azi_2_4 = zeros(FT, cache_dim_ppar+1),
                cache_incl_azi_2_5 = zeros(FT, cache_dim_ppar+1),
                cache_incl_azi_2_6 = zeros(FT, cache_dim_ppar+1),

                cache_incl_azi_1 = zeros(FT, dim_incl, dim_azi),
                cache_incl_azi_2 = zeros(FT, dim_incl, dim_azi),
                cache_incl_azi_3 = zeros(FT, dim_incl, dim_azi),
                cache_incl_azi_4 = zeros(FT, dim_incl, dim_azi),
                cache_incl_azi_5 = zeros(FT, dim_incl, dim_azi),
                cache_incl_azi_6 = zeros(FT, dim_incl, dim_azi),
                cache_incl_azi_7 = zeros(FT, dim_incl, dim_azi),

                solver_nb = NewtonBisectionMethod{FT}(),
                stol_nb   = SolutionTolerance{FT}(eps(FT)*100, 50),
    )
);
