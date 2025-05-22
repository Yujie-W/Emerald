#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the reflectance of the leaf
#
#######################################################################################################################################################################################################
"""

    leaf_ρ(ρ₀::FT, τ₀::FT, ρ₁::FT, τ₁::FT, ρ₂::FT) where {FT}

Return the leaf level reflectance, given
- `ρ₀` reflectance of the first layer with incident radiation
- `τ₀` transmittance of the first layer with incident radiation
- `ρ₁` reflectance of the first layer with isotropic radiation
- `τ₁` transmittance of the first layer with isotropic radiation
- `ρ₂` reflectance of the n-1 layer with isotropic radiation

"""
function leaf_ρ(ρ₀::FT, τ₀::FT, ρ₁::FT, τ₁::FT, ρ₂::FT) where {FT}
    denom = 1 - ρ₁ * ρ₂;

    return ρ₀ + τ₀ * ρ₂ * τ₁ / denom
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute the transmittance of the leaf
#
#######################################################################################################################################################################################################
"""

    leaf_τ(τ₀::FT, ρ₁::FT, ρ₂::FT, τ₂::FT) where {FT}

Return the leaf level reflectance, given
- `τ₀` transmittance of the first layer with incident radiation
- `ρ₁` reflectance of the first layer with isotropic radiation
- `ρ₂` reflectance of the n-1 layer with isotropic radiation
- `τ₂` transmittance of the n-1 layer with isotropic radiation

"""
function leaf_τ(τ₀::FT, ρ₁::FT, ρ₂::FT, τ₂::FT) where {FT}
    denom = 1 - ρ₁ * ρ₂;

    return τ₀ * τ₂ / denom
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to update the leaf reflectance and transmittance within LeafBio
#     2023-Sep-16: save leaf aborption as well
#
#######################################################################################################################################################################################################
"""

    leaf_ρ_τ!(bio::LeafBio{FT}) where {FT}

Update the leaf reflectance and transmittance within `bio`, given
- `bio` LeafBio struct

"""
function leaf_ρ_τ!(bio::LeafBio{FT}) where {FT}
    @. bio.auxil.ρ_leaf = leaf_ρ(bio.auxil.ρ_layer_θ, bio.auxil.τ_layer_θ, bio.auxil.ρ_layer_1, bio.auxil.τ_layer_1, bio.auxil.ρ_layer_2);
    @. bio.auxil.τ_leaf = leaf_τ(bio.auxil.τ_layer_θ, bio.auxil.ρ_layer_1, bio.auxil.ρ_layer_2, bio.auxil.τ_layer_2);
    @. bio.auxil.α_leaf = 1 - bio.auxil.ρ_leaf - bio.auxil.τ_leaf;

    return nothing
end;
