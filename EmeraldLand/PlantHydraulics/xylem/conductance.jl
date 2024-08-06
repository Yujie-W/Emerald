#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Aug-05: add function to return the xylem conductance
#
#######################################################################################################################################################################################################
"""

    xylem_conductance(organ::Union{CanopyLayer{FT}, Leaf{FT}, Root{FT}, Stem{FT}}) where {FT}

Return the xylem conductance, given
- `organ` `CanopyLayer`, `Leaf`, `Root`, or `Stem` type structure

"""
function xylem_conductance end;

xylem_conductance(organ::Union{CanopyLayer{FT}, Leaf{FT}, Root{FT}, Stem{FT}}) where {FT} = xylem_conductance(organ.xylem);

xylem_conductance(xylem::XylemHydraulics{FT}) where {FT} = xylem_conductance(xylem.trait, xylem.state);

xylem_conductance(x_trait::XylemHydraulicsTrait{FT}, x_state::XylemHydraulicsState{FT}) where {FT} = (
    k_max = x_trait.area * x_trait.k_max / x_trait.l;
    r_xyl::FT = 0;

    N = length(x_state.p_history);
    for i in 1:N
        p_mem = x_state.p_history[i];
        k = relative_xylem_k(x_trait.vc, p_mem) * k_max * N;
        r_xyl += 1 / k;
    end;

    return 1 / r_xyl
);
