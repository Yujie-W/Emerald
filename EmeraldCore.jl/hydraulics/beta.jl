#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-30: migrate function from older version StomataModels.jl
#     2022-Jul-01: add method to tune stomatal opening based on relative hydraulic conductance at leaf xylem end
#     2022-Jul-01: add method to tune stomatal opening based on relative hydraulic conductance of the soil
#     2022-Jul-01: add method to tune stomatal opening based on soil potential or leaf pressure
#     2022-Jul-01: fix a typo in function call
#     2022-Jul-12: move function from StomataModels.jl to PlantHydraulics.jl
#     2022-Nov-18: force the beta to be within (0,1]
#
#######################################################################################################################################################################################################
"""

    β_factor(f::Function, vc::AbstractXylemVC{FT}, x_25::FT) where {FT}
    β_factor(f::Function, vc::AbstractSoilVC{FT}, x_25::FT) where {FT}
    β_factor(f::Function, x_25::FT) where {FT}

Return the β factor based on relative conductance or soil potential/pressure, given
- `f` Function to translate relative k to β, for example f(x) = x, f(x) = x², and f(x) = sqrt(x) for x in [0,1]
- `vc` Leaf vulnerability curve or soil vulnerability curve (moisture retention curve)
- `x_25` Leaf xylem pressure corrected to 25 °C, soil water potential corrected to 25 °C (forcing on roots, note that this function may not be useful for plants with salt stress), or soil water
    content.

"""
function β_factor end

β_factor(f::Function, vc::AbstractXylemVC{FT}, x_25::FT) where {FT} = FT(max(eps(FT), min(1, f(relative_hydraulic_conductance(vc, x_25)))));

β_factor(f::Function, vc::AbstractSoilVC{FT}, x_25::FT) where {FT} = FT(max(eps(FT), min(1, f(relative_hydraulic_conductance(vc, true, x_25)))));

β_factor(f::Function, x_25::FT) where {FT} = FT(max(eps(FT), min(1, f(x_25))));

β_factor(sm::AbstractStomataModel{FT}) where {FT} = FT(NaN);

β_factor(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}) where {FT} = sm.β.β₁;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-12: add function to update beta factor for empirical models
#     2022-Jul-12: add methods for MonoElementSPAC
#     2022-Jul-15: rename xylem_flow to flow_in to be more descriptive
#     2022-Oct-19: fix sm.Β to sm.β
#     2022-Oct-20: use add SoilLayer to function variables, because of the removal of SH from RootHydraulics
#     2022-Nov-18: use root K to weigh the beta among root layers
#     2023-Mar-27: weigh the beta among root layers only if flow rate is positive (if all flows are negative, beta = 1)
# Bug fixes
#     2022-Oct-20: fix the issue related to β_factor!(roots, soil, leaves, β, β.PARAM_X) as I forgot to write β_factor! before `()`
#
#######################################################################################################################################################################################################
"""
This function updates the beta factor for SPAC if empirical models are used. The method is meant to support all SPAC defined in EmeraldNamespace.jl:
- `MonoElementSPAC`
- `MultiLayerSPAC`

"""
function β_factor! end

"""

    β_factor!(spac::MonoElementSPAC{FT}) where {FT}

Update the beta factor for the LEAF component in SPAC, given
- `spac` `MonoElementSPAC` type SPAC

"""
β_factor!(spac::MonoElementSPAC{FT}) where {FT} = β_factor!(spac, spac.LEAF.SM);

β_factor!(spac::MonoElementSPAC{FT}, sm::AbstractStomataModel{FT}) where {FT} = nothing;

β_factor!(spac::MonoElementSPAC{FT}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}) where {FT} = β_factor!(spac, sm.β);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}) where {FT} = β_factor!(spac, β, β.PARAM_X);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}, param_x::BetaParameterKleaf) where {FT} = (
    (; LEAF) = spac;

    _f_st = relative_surface_tension(LEAF.t);

    β.β₁ = β_factor(β.FUNC, LEAF.HS.VC, LEAF.HS._p_element[end] / _f_st);

    return nothing
);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}, param_x::BetaParameterKsoil) where {FT} = (
    (; ROOT, SOIL) = spac;

    _f_st = relative_surface_tension(ROOT.t);

    β.β₁ = β_factor(β.FUNC, SOIL.LAYERS[1].VC, ROOT.HS.p_ups / _f_st);

    return nothing
);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}, param_x::BetaParameterPleaf) where {FT} = (
    (; LEAF) = spac;

    _f_st = relative_surface_tension(LEAF.t);

    β.β₁ = β_factor(β.FUNC, LEAF.HS._p_element[end] / _f_st);

    return nothing
);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}, param_x::BetaParameterPsoil) where {FT} = (
    (; ROOT) = spac;

    _f_st = relative_surface_tension(ROOT.t);

    β.β₁ = β_factor(β.FUNC, ROOT.HS.p_ups / _f_st);

    return nothing
);

β_factor!(spac::MonoElementSPAC{FT}, β::BetaFunction{FT}, param_x::BetaParameterΘ) where {FT} = (
    (; ROOT, SOIL) = spac;

    _f_st = relative_surface_tension(ROOT.t);

    β.β₁ = β_factor(β.FUNC, soil_θ(SOIL.LAYERS[1].VC, ROOT.HS.p_ups / _f_st));

    return nothing
);


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

β_factor!(roots::Vector{Root{FT,DIM_XYLEM}}, soil::Soil{FT}, leaves::Leaves2D{FT,DIM_XYLEM}, sm::AbstractStomataModel{FT}) where {FT,DIM_XYLEM} = nothing;

β_factor!(roots::Vector{Root{FT,DIM_XYLEM}}, soil::Soil{FT}, leaves::Leaves2D{FT,DIM_XYLEM}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}) where {FT,DIM_XYLEM} =
    β_factor!(roots, soil, leaves, sm.β);

β_factor!(roots::Vector{Root{FT,DIM_XYLEM}}, soil::Soil{FT}, leaves::Leaves2D{FT,DIM_XYLEM}, β::BetaFunction{FT}) where {FT,DIM_XYLEM} = β_factor!(roots, soil, leaves, β, β.PARAM_X);

β_factor!(roots::Vector{Root{FT,DIM_XYLEM}}, soil::Soil{FT}, leaves::Leaves2D{FT,DIM_XYLEM}, β::BetaFunction{FT}, param_x::BetaParameterKleaf) where {FT,DIM_XYLEM} = (
    _f_st = relative_surface_tension(leaves.t);

    β.β₁ = β_factor(β.FUNC, leaves.HS.VC, leaves.HS._p_element[end] / _f_st);

    return nothing
);

β_factor!(roots::Vector{Root{FT,DIM_XYLEM}}, soil::Soil{FT}, leaves::Leaves2D{FT,DIM_XYLEM}, β::BetaFunction{FT}, param_x::BetaParameterKsoil) where {FT,DIM_XYLEM} = (
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

    if _deno > 0
        β.β₁ = _norm / _deno;
    elseif _sumf < 0
        β.β₁ = 1
    else
        β.β₁ = eps(FT);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT,DIM_XYLEM}}, soil::Soil{FT}, leaves::Leaves2D{FT,DIM_XYLEM}, β::BetaFunction{FT}, param_x::BetaParameterPleaf) where {FT,DIM_XYLEM} = (
    _f_st = relative_surface_tension(leaves.t);

    β.β₁ = β_factor(β.FUNC, leaves.HS._p_element[end] / _f_st);

    return nothing
);

β_factor!(roots::Vector{Root{FT,DIM_XYLEM}}, soil::Soil{FT}, leaves::Leaves2D{FT,DIM_XYLEM}, β::BetaFunction{FT}, param_x::BetaParameterPsoil) where {FT,DIM_XYLEM} = (
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

    if _deno > 0
        β.β₁ = _norm / _deno;
    elseif _sumf < 0
        β.β₁ = 1
    else
        β.β₁ = eps(FT);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT,DIM_XYLEM}}, soil::Soil{FT}, leaves::Leaves2D{FT,DIM_XYLEM}, β::BetaFunction{FT}, param_x::BetaParameterΘ) where {FT,DIM_XYLEM} = (
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

    if _deno > 0
        β.β₁ = _norm / _deno;
    elseif _sumf < 0
        β.β₁ = 1
    else
        β.β₁ = eps(FT);
    end;

    return nothing
);
