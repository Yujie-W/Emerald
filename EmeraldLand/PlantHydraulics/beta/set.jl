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
This function updates the beta factor for SPAC if empirical models are used. The method is meant to support all SPAC defined in EmeraldNamespace.jl:
- `MultiLayerSPAC`

"""
function β_factor! end


"""

    β_factor!(spac::MultiLayerSPAC{FT}) where {FT}

Update the β factor for the LEAVES component in SPAC, given
- `spac` `MultiLayerSPAC` type SPAC

Note that if the β function is based on Kleaf or Pleaf, β factor is taken as that of leaf; if the β function is based on Ksoil, Psoil, or Θ, β is taken as the average weighted by flow rate in each
    root.

"""
β_factor!(spac::MultiLayerSPAC{FT}) where {FT} = (
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
    f_st = relative_surface_tension(leaf.t);

    β.β = β_factor(β.FUNC, leaf.HS.VC, leaf.HS._p_element[end] / f_st);

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, β::BetaFunction{FT}, param_x::BetaParameterKsoil) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    norm = 0;
    sumf = 0;
    denom = 0;
    for i in eachindex(roots)
        f_st = relative_surface_tension(roots[i].auxil.t);
        beta = β_factor(β.FUNC, soils[i].state.vc, roots[i].HS.p_ups / f_st);
        f_in = flow_in(roots[i]);
        kmax = f_in > 0 ? roots[i].HS.AREA * roots[i].HS.K_X / roots[i].HS.L : 0;
        norm += beta * kmax;
        sumf += f_in;
        denom += kmax;
    end;

    @assert !isnan(norm) && !isnan(denom) && !isnan(sumf) "NaN detected in beta calculation!";

    if denom > 0
        β.β = norm / denom;
    elseif sumf < 0
        β.β = 1
    else
        β.β = eps(FT);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, β::BetaFunction{FT}, param_x::BetaParameterPleaf) where {FT} = (
    f_st = relative_surface_tension(leaf.t);

    β.β = β_factor(β.FUNC, leaf.HS._p_element[end] / f_st);

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, β::BetaFunction{FT}, param_x::BetaParameterPsoil) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    norm = 0;
    sumf = 0;
    denom = 0;
    for i in eachindex(roots)
        f_st = relative_surface_tension(roots[i].auxil.t);
        beta = β_factor(β.FUNC, roots[i].HS.p_ups / f_st);
        f_in = flow_in(roots[i]);
        kmax = f_in > 0 ? roots[i].HS.AREA * roots[i].HS.K_X / roots[i].HS.L : 0;
        norm += beta * kmax;
        sumf += f_in;
        denom += kmax;
    end;

    @assert !isnan(norm) && !isnan(denom) && !isnan(sumf) "NaN detected in beta calculation!";

    if denom > 0
        β.β = norm / denom;
    elseif sumf < 0
        β.β = 1
    else
        β.β = eps(FT);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Leaf{FT}, β::BetaFunction{FT}, param_x::BetaParameterΘ) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    norm = 0;
    sumf = 0;
    denom = 0;
    for i in eachindex(roots)
        f_st = relative_surface_tension(roots[i].auxil.t);
        beta = β_factor(β.FUNC, soil_θ(soils[i].state.vc, roots[i].HS.p_ups / f_st));
        f_in = flow_in(roots[i]);
        kmax = f_in > 0 ? roots[i].HS.AREA * roots[i].HS.K_X / roots[i].HS.L : 0;
        norm += beta * kmax;
        sumf += f_in;
        denom += kmax;
    end;

    @assert !isnan(norm) && !isnan(denom) && !isnan(sumf) "NaN detected in beta calculation!";

    if denom > 0
        β.β = norm / denom;
    elseif sumf < 0
        β.β = 1
    else
        β.β = eps(FT);
    end;

    return nothing
);
