# This file contains function to set beta factor for the stomatal models

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-12: add function to update beta factor for empirical models
#     2022-Jul-15: rename xylem_flow to flow_in to be more descriptive
#     2022-Oct-19: fix sm.Β to sm.β
#     2022-Oct-20: use add SoilLayer to function variables, because of the removal of SH from RootHydraulics
#     2022-Nov-18: use root K to weigh the beta among root layers
#     2023-Mar-27: weigh the beta among root layers only if flow rate is positive (if all flows are negative, beta = 1)
#     2023-Aug-27: add nan check for beta calculation of empirical models
# Bug fixes
#     2022-Oct-20: fix the issue related to β_factor!(roots, soil, leaf, β, β.PARAM_X) as I forgot to write β_factor! before `()`
#
#######################################################################################################################################################################################################
"""

    β_factor!(spac::BulkSPAC{FT}) where {FT}

Update the β factor for the LEAVES component in SPAC, given
- `spac` `BulkSPAC` type SPAC

Note that if the β function is based on Kleaf or Pleaf, β factor is taken as that of leaf; if the β function is based on Ksoil, Psoil, or Θ, β is taken as the average weighted by flow rate in each
    root.

"""
function β_factor! end;

β_factor!(spac::BulkSPAC{FT}) where {FT} = (
    (; LEAVES, ROOTS, SOILS) = spac;

    for i in eachindex(LEAVES)
        β_factor!(ROOTS, SOILS, LEAVES[i], LEAVES[i].flux.state.stomatal_model);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, sm::AbstractStomataModel{FT}) where {FT} = nothing;

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}) where {FT} =
    β_factor!(roots, soils, leaf, sm.β);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, β::BetaFunction{FT}) where {FT} = β_factor!(roots, soils, leaf, β, β.PARAM_X);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, β::BetaFunction{FT}, param_x::BetaParameterKleaf) where {FT} = (
    f_st = relative_surface_tension(leaf.energy.auxil.t);

    leaf.flux.auxil.β = β_factor(β.FUNC, relative_xylem_k(leaf.xylem.state.vc, leaf.xylem.auxil.pressure[end] / f_st));

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, β::BetaFunction{FT}, param_x::BetaParameterKsoil) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    norm = 0;
    sumf = 0;
    denom = 0;
    for i in eachindex(roots)
        beta = β_factor(β.FUNC, soils[i].auxil.k);
        f_in = flow_in(roots[i]);
        kmax = f_in > 0 ? roots[i].xylem.state.area * roots[i].xylem.state.k_max / roots[i].xylem.state.l : 0;
        norm += beta * kmax;
        sumf += f_in;
        denom += kmax;
    end;

    @assert !isnan(norm) && !isnan(denom) && !isnan(sumf) "NaN detected in beta calculation!";

    if denom > 0
        leaf.flux.auxil.β = norm / denom;
    elseif sumf < 0
        leaf.flux.auxil.β = 1
    else
        leaf.flux.auxil.β = eps(FT);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, β::BetaFunction{FT}, param_x::BetaParameterPleaf) where {FT} = (
    leaf.flux.auxil.β = β_factor(β.FUNC, leaf.xylem.auxil.pressure[end]);

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, β::BetaFunction{FT}, param_x::BetaParameterPsoil) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    norm = 0;
    sumf = 0;
    denom = 0;
    for i in eachindex(roots)
        beta = β_factor(β.FUNC, soils[i].auxil.ψ);
        f_in = flow_in(roots[i]);
        kmax = f_in > 0 ? roots[i].xylem.state.area * roots[i].xylem.state.k_max / roots[i].xylem.state.l : 0;
        norm += beta * kmax;
        sumf += f_in;
        denom += kmax;
    end;

    @assert !isnan(norm) && !isnan(denom) && !isnan(sumf) "NaN detected in beta calculation!";

    if denom > 0
        leaf.flux.auxil.β = norm / denom;
    elseif sumf < 0
        leaf.flux.auxil.β = 1
    else
        leaf.flux.auxil.β = eps(FT);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, β::BetaFunction{FT}, param_x::BetaParameterΘ) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    norm = 0;
    sumf = 0;
    denom = 0;
    for i in eachindex(roots)
        beta = β_factor(β.FUNC, soils[i].state.θ);
        f_in = flow_in(roots[i]);
        kmax = f_in > 0 ? roots[i].xylem.state.area * roots[i].xylem.state.k_max / roots[i].xylem.state.l : 0;
        norm += beta * kmax;
        sumf += f_in;
        denom += kmax;
    end;

    @assert !isnan(norm) && !isnan(denom) && !isnan(sumf) "NaN detected in beta calculation!";

    if denom > 0
        leaf.flux.auxil.β = norm / denom;
    elseif sumf < 0
        leaf.flux.auxil.β = 1
    else
        leaf.flux.auxil.β = eps(FT);
    end;

    return nothing
);
