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
#     2022-Oct-20: fix the issue related to β_factor!(roots, soil, leaves, β, β.PARAM_X) as I forgot to write β_factor! before `()`
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
    (; LEAVES, ROOTS, SOIL) = spac;

    for _i in eachindex(LEAVES)
        β_factor!(ROOTS, SOIL, LEAVES[_i], LEAVES[_i].SM);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soil::Soil{FT}, leaves::Leaves2D{FT}, sm::AbstractStomataModel{FT}) where {FT} = nothing;

β_factor!(roots::Vector{Root{FT}}, soil::Soil{FT}, leaves::Leaves2D{FT}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}) where {FT} =
    β_factor!(roots, soil, leaves, sm.β);

β_factor!(roots::Vector{Root{FT}}, soil::Soil{FT}, leaves::Leaves2D{FT}, β::BetaFunction{FT}) where {FT} = β_factor!(roots, soil, leaves, β, β.PARAM_X);

β_factor!(roots::Vector{Root{FT}}, soil::Soil{FT}, leaves::Leaves2D{FT}, β::BetaFunction{FT}, param_x::BetaParameterKleaf) where {FT} = (
    _f_st = relative_surface_tension(leaves.t);

    β.β = β_factor(β.FUNC, leaves.HS.VC, leaves.HS._p_element[end] / _f_st);

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soil::Soil{FT}, leaves::Leaves2D{FT}, β::BetaFunction{FT}, param_x::BetaParameterKsoil) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    _norm = 0;
    _deno = 0;
    _sumf = 0;
    for _i in eachindex(roots)
        _f_st = relative_surface_tension(roots[_i].t);
        _beta = β_factor(β.FUNC, soil.LAYERS[_i].VC, roots[_i].HS.p_ups / _f_st);
        _f_in = flow_in(roots[_i]);
        _kmax = _f_in > 0 ? roots[_i].HS.AREA * roots[_i].HS.K_X / roots[_i].HS.L : 0;
        _norm += _beta * _kmax;
        _deno += _kmax;
        _sumf += _f_in;
    end;

    @assert !isnan(_norm) && !isnan(_deno) && !isnan(_sumf) "NaN detected in beta calculation!";

    if _deno > 0
        β.β = _norm / _deno;
    elseif _sumf < 0
        β.β = 1
    else
        β.β = eps(FT);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soil::Soil{FT}, leaves::Leaves2D{FT}, β::BetaFunction{FT}, param_x::BetaParameterPleaf) where {FT} = (
    _f_st = relative_surface_tension(leaves.t);

    β.β = β_factor(β.FUNC, leaves.HS._p_element[end] / _f_st);

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soil::Soil{FT}, leaves::Leaves2D{FT}, β::BetaFunction{FT}, param_x::BetaParameterPsoil) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    _norm = 0;
    _deno = 0;
    _sumf = 0;
    for _i in eachindex(roots)
        _f_st = relative_surface_tension(roots[_i].t);
        _beta = β_factor(β.FUNC, roots[_i].HS.p_ups / _f_st);
        _f_in = flow_in(roots[_i]);
        _kmax = _f_in > 0 ? roots[_i].HS.AREA * roots[_i].HS.K_X / roots[_i].HS.L : 0;
        _norm += _beta * _kmax;
        _deno += _kmax;
        _sumf += _f_in;
    end;

    @assert !isnan(_norm) && !isnan(_deno) && !isnan(_sumf) "NaN detected in beta calculation!";

    if _deno > 0
        β.β = _norm / _deno;
    elseif _sumf < 0
        β.β = 1
    else
        β.β = eps(FT);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soil::Soil{FT}, leaves::Leaves2D{FT}, β::BetaFunction{FT}, param_x::BetaParameterΘ) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    _norm = 0;
    _deno = 0;
    _sumf = 0;
    for _i in eachindex(roots)
        _f_st = relative_surface_tension(roots[_i].t);
        _beta = β_factor(β.FUNC, soil_θ(soil.LAYERS[_i].VC, roots[_i].HS.p_ups / _f_st));
        _f_in = flow_in(roots[_i]);
        _kmax = _f_in > 0 ? roots[_i].HS.AREA * roots[_i].HS.K_X / roots[_i].HS.L : 0;
        _norm += _beta * _kmax;
        _deno += _kmax;
        _sumf += _f_in;
    end;

    @assert !isnan(_norm) && !isnan(_deno) && !isnan(_sumf) "NaN detected in beta calculation!";

    if _deno > 0
        β.β = _norm / _deno;
    elseif _sumf < 0
        β.β = 1
    else
        β.β = eps(FT);
    end;

    return nothing
);
