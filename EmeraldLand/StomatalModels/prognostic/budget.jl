# This file contains functions to update prognostic stomatal conductance

#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add new method to update stomatal conductance prognostically
#     2022-Jul-12: add method to update g for SPAC
#     2023-Mar-11: do nothing if LAI == 0
#     2023-Sep-14: add root disconnection control
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Update stomatal conductance for H₂O based on computed ∂g∂t, given
- `spac` `MultiLayerSPAC` type struct
- `δt` Time step length `[s]`

"""
stomatal_conductance!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT} = (
    # if lai = 0 or roots are not connected, do nothing
    if spac.CANOPY.structure.state.lai == 0 || !spac._root_connection
        return nothing
    end;

    (; LEAVES) = spac;

    for leaf in LEAVES
        stomatal_conductance!(leaf, δt);
    end;

    return nothing
);

stomatal_conductance!(leaf::Leaf{FT}, δt::FT) where {FT} = (
    leaf.flux.state.g_H₂O_s_shaded += leaf.flux.auxil.∂g∂t_shaded * δt;
    for i in eachindex(leaf.flux.state.g_H₂O_s_sunlit)
        leaf.flux.state.g_H₂O_s_sunlit[i] += leaf.flux.auxil.∂g∂t_sunlit[i] * δt;
    end;

    return nothing
);
