"""
Dimensions and angles of the SPAC model
"""
Base.@kwdef mutable struct SPACDimensions{FT<:AbstractFloat}
    # dimensions
    "Dimension of azimuth angles"
    DIM_AZI::Int = 36
    "Dimension of inclination angles"
    DIM_INCL::Int = 9
    "Number of sunlit PPAR bins for all the layers (to speed up the computation; 0 for one leaf model)"
    DIM_PPAR_BINS::Union{Int,Nothing} = nothing
    "Dimension of xylem slices of leaf, stem, and root; xylem capaciatance of stem and root"
    DIM_XYLEM::Int = 5

    # angles
    "Mean azimuth angles `[°]`"
    Θ_AZI::Vector{FT} = collect(FT, range(0, 360; length=DIM_AZI+1))[1:end-1] .+ 360 / DIM_AZI / 2
    "Bounds of inclination angles `[°]`"
    Θ_INCL_BNDS::Matrix{FT} = FT[ collect(FT, range(0, 90; length=DIM_INCL+1))[1:end-1] collect(FT, range(0, 90; length=DIM_INCL+1))[2:end] ]
    "Mean inclination angles `[°]`"
    Θ_INCL::Vector{FT} = FT[ (Θ_INCL_BNDS[i,1] + Θ_INCL_BNDS[i,2]) / 2 for i in 1:DIM_INCL ]
end;
