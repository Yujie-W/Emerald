"""

    prescribe_ps_traits!(leaf::Leaf{FT}; args...) where {FT}
    prescribe_ps_traits!(spac::BulkSPAC{FT}; args...) where {FT}

Prescribe the photosynthetic traits for a single leaf, given
- `leaf` `Leaf` or `Leaf` object, keyword arguments:
    - `b6f` b6f content. Optional, default is nothing
    - `jmax` Jmax25. Optional, default is nothing
    - `rd` Dark respiration rate. Optional, default is nothing
    - `vcmax` Vcmax25. Optional, default is nothing
    - `vpmax` Vpmax25. Optional, default is nothing
- `spac` BulkSPAC object, keyword arguments:
    - `vertical_expo` Exponential tuning factor to adjust Vcmax25. Optional, default is nothing

"""
function prescribe_ps_traits! end;

# Prescribe the variables for single leaf, suggest to do this only for top canopy
prescribe_ps_traits!(leaf::Leaf{FT}; args...) where {FT} = prescribe_ps_traits!(leaf.photosystem.trait; args...);

prescribe_ps_traits!(
            pst::GeneralC3Trait{FT};
            b6f::Union{Nothing,Number} = nothing,
            jmax::Union{Nothing,Number} = nothing,
            rd::Union{Nothing,Number} = nothing,
            vcmax::Union{Nothing,Number} = nothing,
            vpmax::Union{Nothing,Number} = nothing) where {FT} = (
    isnothing(vcmax) ? nothing : pst.v_cmax25 = vcmax;
    isnothing(jmax)  ? nothing : pst.j_max25  = jmax;
    isnothing(b6f)   ? nothing : pst.b₆f      = b6f;
    isnothing(rd)    ? nothing : pst.r_d25    = rd;

    return nothing
);

prescribe_ps_traits!(
            pst::GeneralC4Trait{FT};
            b6f::Union{Nothing,Number} = nothing,
            jmax::Union{Nothing,Number} = nothing,
            rd::Union{Nothing,Number} = nothing,
            vcmax::Union{Nothing,Number} = nothing,
            vpmax::Union{Nothing,Number} = nothing) where {FT} = (
    isnothing(vcmax) ? nothing : pst.v_cmax25 = vcmax;
    isnothing(vpmax) ? nothing : pst.v_pmax25 = vpmax;
    isnothing(rd)    ? nothing : pst.r_d25    = rd;

    return nothing
);

# Method to apply the exponential tuning factor to Vcmax25...
prescribe_ps_traits!(spac::BulkSPAC{FT}; args...) where {FT} = prescribe_ps_traits!(spac, spac.plant.leaves[end].photosystem.trait; args...);

prescribe_ps_traits!(spac::BulkSPAC{FT}, ::GeneralC3Trait{FT}; vertical_expo::Union{Nothing,Number} = nothing) where {FT} = (
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # update vertical profiles
    for irt in 1:n_layer
        ilf = n_layer - irt + 1;
        ratio = isnothing(vertical_expo) ? 1 : exp(-vertical_expo * sum(view(can_str.trait.δlai,1:irt-1)));
        leaf = leaves[ilf];
        leaf.photosystem.trait.v_cmax25 = leaves[end].photosystem.trait.v_cmax25 * ratio;
        leaf.photosystem.trait.j_max25 = leaves[end].photosystem.trait.j_max25 * ratio;
        leaf.photosystem.trait.r_d25 = leaves[end].photosystem.trait.r_d25 * ratio;
        leaf.photosystem.trait.b₆f = leaves[end].photosystem.trait.b₆f * ratio;
    end;

    return nothing
);

prescribe_ps_traits!(spac::BulkSPAC{FT}, ::GeneralC4Trait{FT}; vertical_expo::Union{Nothing,Number} = nothing) where {FT} = (
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # update vertical profiles
    for irt in 1:n_layer
        ilf = n_layer - irt + 1;
        ratio = isnothing(vertical_expo) ? 1 : exp(-vertical_expo * sum(view(can_str.trait.δlai,1:irt-1)));
        leaf = leaves[ilf];
        leaf.photosystem.trait.v_cmax25 = leaves[end].photosystem.trait.v_cmax25 * ratio;
        leaf.photosystem.trait.v_pmax25 = leaves[end].photosystem.trait.v_pmax25 * ratio;
        leaf.photosystem.trait.r_d25 = leaves[end].photosystem.trait.r_d25 * ratio;
    end;

    return nothing
);
