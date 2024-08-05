#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Aug-05: add function to return the xylem conductance
#
#######################################################################################################################################################################################################
"""

    xylem_conductance(organ::Union{Leaf{FT}, Root{FT}, Stem{FT}}) where {FT}
    xylem_conductance(xylem::XylemHydraulics{FT}) where {FT}

Return the xylem conductance, given
- `organ` `Leaf` `Root` `Stem` type struct
- `xylem` `XylemHydraulics` type struct

"""
function xylem_conductance end;

xylem_conductance(layer::CanopyLayer{FT}) where {FT} = xylem_conductance(layer.xylem);

xylem_conductance(leaf::Leaf{FT}) where {FT} = xylem_conductance(leaf.xylem);

xylem_conductance(root::Root{FT}) where {FT} = xylem_conductance(root.xylem);

xylem_conductance(stem::Stem{FT}) where {FT} = xylem_conductance(stem.xylem);

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
