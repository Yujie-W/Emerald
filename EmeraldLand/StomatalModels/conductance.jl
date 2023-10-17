




#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jul-07: add new function
#     2022-Jul-11: deflate documentations
#
#######################################################################################################################################################################################################
"""
This function updates stomatal conductance for H₂O and CO₂. Supported functionalities are
- Update conductance for H₂O prognostically
- Update conductance for CO₂ based on that for H₂O

"""
function stomatal_conductance! end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-12: add new method to compute marginal stomatal conductance increase
#     2022-Jul-12: add new method to update marginal stomatal conductance for SPAC
#     2023-Sep-14: add root disconnection control
# To do
#     TODO: be careful with the β here (need to used the value stored in empirical stomtal model)
#     TODO: use ∂gₙ∂t for nighttime conditions
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance!(spac::MultiLayerSPAC{FT}; β::FT = FT(1)) where {FT}

Update marginal stomatal conductance, given
- `spac` `MultiLayerSPAC` type struct
- `β` Tuning factor

"""
stomatal_conductance!(spac::MultiLayerSPAC{FT}; β::FT = FT(1)) where {FT} = (
    (; AIRS, CANOPY, LEAVES, LEAVES_INDEX) = spac;

    # if lai = 0 or roots are not connected, do nothing
    # TODO: redo this later to foce dgdt to 0
    if CANOPY.structure.state.lai == 0 || !spac._root_connection
        return nothing
    end;

    for i in eachindex(LEAVES_INDEX)
        stomatal_conductance!(LEAVES[i], AIRS[LEAVES_INDEX[i]]; β = β);
    end;

    return nothing
);

stomatal_conductance!(leaf::Leaf{FT}, air::AirLayer{FT}; β::FT = FT(1)) where {FT} = (
    leaf.flux.auxil.∂g∂t_shaded = ∂g∂t(leaf, air; β = β);
    for i in eachindex(leaf.flux.auxil.∂g∂t_sunlit)
        leaf.flux.auxil.∂g∂t_sunlit[i] = ∂g∂t(leaf, air, i; β = β);
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add new method to update stomatal conductance for CO₂ based on that of H₂O
#     2022-Jul-07: add new method to update stomatal conductance prognostically
#     2022-Jul-12: move ∂g∂t to another method
#     2022-Jul-12: add method to update g for SPAC
#     2022-Jul-26: limit g in range after updating stomatal conductance
#     2023-Mar-11: do nothing if LAI == 0
#     2023-Mar-13: move some methods as stomatal_conductance_profile!
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
    (; CANOPY, LEAVES) = spac;

    # if lai = 0 or roots are not connected, do nothing
    if CANOPY.structure.state.lai == 0 || !spac._root_connection
        return nothing
    end;

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
    limit_stomatal_conductance!(leaf);
    stomatal_conductance_profile!(leaf);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2023-Mar-13: add function to update stomatal conductance profile based on gs and gb
#
#######################################################################################################################################################################################################
"""

    stomatal_conductance_profile!(spac::MultiLayerSPAC{FT}) where {FT}

Update stomatal conductance for CO₂ based on that for H₂O, given
- `spac` `MultiLayerSPAC` type struct

"""
function stomatal_conductance_profile! end;

stomatal_conductance_profile!(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; CANOPY, LEAVES) = spac;

    if CANOPY.structure.state.lai == 0
        return nothing
    end;

    for leaf in LEAVES
        stomatal_conductance_profile!(leaf);
    end;

    return nothing
);

stomatal_conductance_profile!(leaf::Leaf{FT}) where {FT} = (
    leaf.flux.auxil.g_CO₂_shaded = 1 / (1 / leaf.flux.auxil.g_CO₂_b + FT(1.6) / leaf.flux.state.g_H₂O_s_shaded);
    for i in eachindex(leaf.flux.state.g_H₂O_s_sunlit)
        leaf.flux.auxil.g_CO₂_sunlit[i] = 1 / (1 / leaf.flux.auxil.g_CO₂_b + FT(1.6) / leaf.flux.state.g_H₂O_s_sunlit[i]);
    end;

    return nothing
);
