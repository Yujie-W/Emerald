# This file contains functions to construct soil hydraulic properties from other soil hydraulic properties or data

#######################################################################################################################################################################################################
#
# Changes made to this constructor
# General
#     2021-Sep-30: move this function out of BrooksCorey struct as an external method for the constructor (avoid dependency on ConstrainedRootSolvers)
#     2022-Apr-19: fix documentation
#     2022-Jul-15: BrooksCorey field changed, modify the constructor accordingly
#
#######################################################################################################################################################################################################
"""

    BrooksCorey{FT}(vg::VanGenuchten{FT}) where {FT}

A constructor for BrooksCorey to create BrooksCorey type soil from VanGenuchten type, given
- `vg` `VanGenuchten` type soil water retention curve
"""
BrooksCorey{FT}(vg::VanGenuchten{FT}) where {FT} = (
    bc = BrooksCorey{FT}(K_MAX = vg.K_MAX, B = 1, TYPE = vg.TYPE, Ψ_SAT = 0.001, Θ_SAT = vg.Θ_SAT, Θ_RES = vg.Θ_RES);

    # generate data to fit
    Θs   = range(vg.Θ_RES+FT(1e-2); stop=vg.Θ_SAT-FT(1e-2), length=30);
    ψ_vG = -1 .* soil_ψ_25.([vg], Θs);
    ψ_BC = similar(ψ_vG);
    ψ_DF = similar(ψ_vG);

    # function to fit BrooksCorey parameters
    @inline fit_bc(x) = (
        bc.B = x[1];
        bc.Ψ_SAT = x[2];
        ψ_BC .= -1 .* soil_ψ_25.([bc], Θs);
        ψ_DF .= (log.(ψ_BC) .- log.(ψ_vG)) .^ 2;

        return -sum(ψ_DF)
    );

    st  = SolutionToleranceND{FT}([1e-3, 1e-6], 30);
    ms  = ReduceStepMethodND{FT}(x_mins = FT[1e-3, 1e-6], x_maxs = FT[ 100, 1000], x_inis = [(2*vg.N-1) / (vg.N-1), 1 / (vg.α)], Δ_inis = FT[0.1, 1e-3]);
    sol = find_peak(fit_bc, ms, st);
    bc.B = sol[1];
    bc.Ψ_SAT = sol[2];

    return bc
);
