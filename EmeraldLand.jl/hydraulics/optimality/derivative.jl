#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jul-08: migrate function from StomataModels.jl
#     2022-Jul-08: add method for LeafHydraulics
#     2022-Jul-08: add methods for Leaves2D
#     2022-Jul-08: add option δe to be more general
#     2023-Aug-27: add nan check
#
#######################################################################################################################################################################################################
"""

    ∂E∂P(lf::Leaves2D{FT}, flow::FT; δe::FT = FT(1e-7)) where {FT}

Return the marginal hydraulic conductance, given
- `lf` `Leaves2D` type struct
- `flow` Flow rate through the leaf xylem `[mol m⁻² s⁻¹]`
- `δe` Incremental flow rate, default is 1e-7

"""
function ∂E∂P end

∂E∂P(lf::Leaves2D{FT}, flow::FT; δe::FT = FT(1e-7)) where {FT} = ∂E∂P(lf.HS, flow, lf.t; δe = δe);

∂E∂P(xylem::XylemHydraulics{FT}, flow::FT, t::FT; δe::FT = FT(1e-7)) where {FT} = (
    @assert δe != 0 "δe must not be 0";

    p1 = xylem_end_pressure(xylem, flow, t);
    p2 = xylem_end_pressure(xylem, flow + δe, t);
    dedp = -δe / (p2 - p1);

    if isnan(dedp)
        @info "Debugging" p1 p2 flow t δe;
        error("NaN detected when computing ∂E∂P!");
    end;

    return dedp
);
