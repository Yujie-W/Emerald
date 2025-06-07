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
#     2024-Feb-29: add LAI = 0 controller
#     2024-Sep-03: use state.asap to check the xylem status (<= 0 means the xylem is dead)
#     2025-Jun-05: when computing beta from soil water content or potential, use eps when there is ice
#     2025-Jun-07: set beta to eps(FT) when there is ice in the root zone
# Bug fixes
#     2022-Oct-20: fix the issue related to β_factor!(roots, soil, leaf, β, β.PARAM_X) as I forgot to write β_factor! before `()`
# To-dos
#     TODO: add ice control for the cases when plant hydraulics takes effect
#
#######################################################################################################################################################################################################
"""

    β_factor!(spac::BulkSPAC{FT}) where {FT}

Update the β factor for the leaves in SPAC, given
- `spac` `BulkSPAC` type SPAC

Note that if the β function is based on Kleaf or Pleaf, β factor is taken as that of leaf; if the β function is based on Ksoil, Psoil, or Θ, β is taken as the average weighted by flow rate in each
    root.

"""
function β_factor! end;

β_factor!(spac::BulkSPAC{FT}) where {FT} = (
    if spac.canopy.structure.trait.lai <= 0
        return nothing
    end;

    # update the beta factor for the leaves only if the LAI is positive
    leaves = spac.plant.leaves;
    roots = spac.plant.roots;
    soils = spac.soils;

    for i in eachindex(leaves)
        β_factor!(roots, soils, leaves[i], leaves[i].flux.trait.stomatal_model);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, sm::AbstractStomataModel{FT}) where {FT} = nothing;

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}) where {FT} = (
    if leaf.xylem.state.asap <= 0
        return nothing
    end;

    # update the beta factor for the leaf only if the LAI is positive
    β_factor!(roots, soils, leaf, sm.β);

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, β::BetaFunction{FT}) where {FT} = β_factor!(roots, soils, leaf, β, β.PARAM_X);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, β::BetaFunction{FT}, param_x::BetaParameterKleaf) where {FT} = (
    f_st = relative_surface_tension(leaf.energy.s_aux.t);

    leaf.flux.auxil.β = β_factor(β.FUNC, relative_xylem_k(leaf.xylem.trait.vc, leaf.xylem.auxil.pressure[end] / f_st));

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, β::BetaFunction{FT}, param_x::BetaParameterKsoil) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    norm = 0;
    sumf = 0;
    denom = 0;
    betas = false;
    for i in eachindex(roots)
        beta = β_factor(β.FUNC, soils[i].s_aux.k);
        f_in = flow_in(roots[i]);
        kmax = f_in > 0 ? roots[i].xylem.state.asap * roots[i].xylem.trait.k_max / roots[i].xylem.trait.l : 0;
        norm += beta * kmax;
        sumf += f_in;
        denom += kmax;

        # if ice > 0.05, then we need to use eps(FT) for the beta
        if soils[i].state.θ_ice > 0.05
            betas = true;
        end;
    end;

    @assert !isnan(norm) && !isnan(denom) && !isnan(sumf) "NaN detected in beta calculation!";

    if denom > 0
        leaf.flux.auxil.β = norm / denom;
    elseif sumf < 0
        leaf.flux.auxil.β = 1;
    else
        leaf.flux.auxil.β = eps(FT);
    end;

    # if ice exists, then we need to use eps(FT) for the beta
    if betas
        leaf.flux.auxil.β = eps(FT);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, β::BetaFunction{FT}, param_x::BetaParameterPleaf) where {FT} = (
    leaf.flux.auxil.β = β_factor(β.FUNC, leaf.xylem.auxil.pressure[end]);

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, β::BetaFunction{FT}, param_x::BetaParameterPsoil) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    norm = 0;
    sumf = 0;
    denom = 0;
    betas = false;
    for i in eachindex(roots)
        beta = β_factor(β.FUNC, soils[i].s_aux.ψ);
        f_in = flow_in(roots[i]);
        kmax = f_in > 0 ? roots[i].xylem.state.asap * roots[i].xylem.trait.k_max / roots[i].xylem.trait.l : 0;
        norm += beta * kmax;
        sumf += f_in;
        denom += kmax;

        # if ice > 0.05, then we need to use eps(FT) for the beta
        if soils[i].state.θ_ice > 0.05
            betas = true;
        end;
    end;

    @assert !isnan(norm) && !isnan(denom) && !isnan(sumf) "NaN detected in beta calculation!";

    if denom > 0
        leaf.flux.auxil.β = norm / denom;
    elseif sumf < 0
        leaf.flux.auxil.β = 1
    else
        leaf.flux.auxil.β = eps(FT);
    end;

    # if ice exists, then we need to use eps(FT) for the beta
    if betas
        leaf.flux.auxil.β = eps(FT);
    end;

    return nothing
);

β_factor!(roots::Vector{Root{FT}}, soils::Vector{SoilLayer{FT}}, leaf::Union{CanopyLayer{FT}, Leaf{FT}}, β::BetaFunction{FT}, param_x::BetaParameterΘ) where {FT} = (
    # weigh the beta by root Kmax for the roots with positive flow
    norm = 0;
    sumf = 0;
    denom = 0;
    betas = false;
    for i in eachindex(roots)
        beta = β_factor(β.FUNC, soils[i].state.θ);
        f_in = flow_in(roots[i]);
        kmax = f_in > 0 ? roots[i].xylem.state.asap * roots[i].xylem.trait.k_max / roots[i].xylem.trait.l : 0;
        norm += beta * kmax;
        sumf += f_in;
        denom += kmax;

        # if ice > 0.05, then we need to use eps(FT) for the beta
        if soils[i].state.θ_ice > 0.05
            betas = true;
        end;
    end;

    @assert !isnan(norm) && !isnan(denom) && !isnan(sumf) "NaN detected in beta calculation!";

    if denom > 0
        leaf.flux.auxil.β = norm / denom;
    elseif sumf < 0
        leaf.flux.auxil.β = 1
    else
        leaf.flux.auxil.β = eps(FT);
    end;

    # if ice exists, then we need to use eps(FT) for the beta
    if betas
        leaf.flux.auxil.β = eps(FT);
    end;

    return nothing
);
