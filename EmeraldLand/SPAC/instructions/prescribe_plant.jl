# This file contains functions to prescribe plant parameters

#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Oct-19: add method to update or prescribe cab, car, lai, Vcmax and Jmax TD, t_leaf, vcmax profile
#     2022-Nov-21: fix a bug related to Vcmax profile (no global simulations are impacted)
#     2023-May-11: add ci to the option list
#     2023-May-19: use δlai per canopy layer
#     2023-Aug-25: add option to set up hydraulic conductance profiles for root, trunk, branches, and leaves
#     2023-Aug-27: fix a typo in the computation of k profiles (reverse the denominator and numerator)
#     2023-Oct-02: run energy initialization when LAI or t_leaf is updated
#     2023-Oct-18: recalculate canopy structural parameters when LAI, cab, car, ci is updated
#     2024-Jan-11: add option to prescribe sai
#     2024-Feb-08: fix the issue in Vcmax profiles (was opposite)
#     2024-Feb-08: add support to C3State
#     2024-Feb-27: run dull_aux! when any of the canopy structural parameters or leaf pigments is updated
#     2024-Feb-28: set minimum LAI to 0 if a negative value is prescribed
#     2024-Jul-24: add leaf shedded flag
#
#######################################################################################################################################################################################################
"""

    prescribe_traits!(
                config::SPACConfiguration{FT},
                spac::BulkSPAC{FT};
                b6f::Union{Number,Nothing} = nothing,
                cab::Union{Number,Nothing} = nothing,
                car::Union{Number,Nothing} = nothing,
                ci::Union{Number,Nothing} = nothing,
                jmax::Union{Number,Nothing} = nothing,
                kmax::Union{Number,Tuple,Nothing} = nothing,
                lai::Union{Number,Nothing} = nothing,
                rd::Union{Number,Nothing} = nothing,
                sai::Union{Number,Nothing} = nothing,
                t_clm::Union{Number,Nothing} = nothing,
                t_leaf::Union{Number,Nothing} = nothing,
                vcmax::Union{Number,Nothing} = nothing,
                vertical_expo::Union{Number,Nothing} = nothing,
                vpmax::Union{Number,Nothing} = nothing) where {FT}

Update the physiological parameters of the SPAC, given
- `config` Configuration for `BulkSPAC`
- `spac` Soil plant air continuum
- `b6f` b6f content for C3Cyto model at the top of canopy. Optional, default is nothing
- `cab` Chlorophyll content. Optional, default is nothing
- `car` Carotenoid content. Optional, default is nothing
- `ci` Clumping index. Optional, default is nothing
- `jmax` Jmax25 at the top of canopy. Optional, default is nothing
- `kmax` Maximum hydraulic conductance. Optional, default is nothing
- `lai` Leaf area index. Optional, default is nothing
- `rd` Dark respiration rate at the top of canopy. Optional, default is nothing
- `sai` Stem area index. Optional, default is nothing
- `t_clm` Moving average temperature to update Vcmax and Jmax temperature dependencies. Optional, default is nothing
- `t_leaf` Leaf temperature. Optional, default is nothing
- `vcmax` Vcmax25 at the top of canopy. Optional, default is nothing
- `vertical_expo` Exponential tuning factor to adjust Vcmax25. Optional, default is nothing
- `vpmax` Vpmax25 at the top of canopy. Optional, default is nothing

"""
function prescribe_traits!(
            config::SPACConfiguration{FT},
            spac::BulkSPAC{FT};
            b6f::Union{Number,Nothing} = nothing,
            cab::Union{Number,Nothing} = nothing,
            car::Union{Number,Nothing} = nothing,
            ci::Union{Number,Nothing} = nothing,
            jmax::Union{Number,Nothing} = nothing,
            kmax::Union{Number,Tuple,Nothing} = nothing,
            lai::Union{Number,Nothing} = nothing,
            rd::Union{Number,Nothing} = nothing,
            sai::Union{Number,Nothing} = nothing,
            t_clm::Union{Number,Nothing} = nothing,
            t_leaf::Union{Number,Nothing} = nothing,
            vcmax::Union{Number,Nothing} = nothing,
            vertical_expo::Union{Number,Nothing} = nothing,
            vpmax::Union{Number,Nothing} = nothing,
) where {FT}
    branches = spac.plant.branches;
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    roots = spac.plant.roots;
    sbulk = spac.soil_bulk;
    trunk = spac.plant.trunk;
    n_layer = length(leaves);

    #
    # plant traits
    #
    # update chlorophyll and carotenoid contents (if leaf shedding flag is not true)
    if !isnothing(cab) && !spac.plant._leaf_shedded
        for leaf in leaves
            leaf.bio.trait.cab = cab;
        end;
    end;

    if !isnothing(car) && !spac.plant._leaf_shedded
        for leaf in leaves
            leaf.bio.trait.car = car;
        end;
    end;

    # update kmax (if leaf shedding flag is not true)
    if !isnothing(kmax) && !spac.plant._leaf_shedded
        # set up the kmax assuming 50% resistance in root, 25% in stem, and 25% in leaves
        ks = if kmax isa Number
            trunk_percent = trunk.xylem.trait.Δh / (trunk.xylem.trait.Δh + branches[end].xylem.trait.Δh);
            (2 * kmax, 4 * kmax / trunk_percent, 4 * kmax / (1 - trunk_percent), 4 * kmax)
        else
            @assert length(kmax) == 4 "kmax must be a number or a tuple of length 4";
            kmax
        end;

        # partition kmax into the roots based on xylem area
        for root in roots
            # root.xylem.trait.k_max = root.xylem.trait.area / trunk.xylem.trait.area * ks[1] * root.xylem.trait.l / root.xylem.trait.area;
            root.xylem.trait.k_max = ks[1] * root.xylem.trait.l / trunk.xylem.trait.area;
        end;
        trunk.xylem.trait.k_max = ks[2] * trunk.xylem.trait.l / trunk.xylem.trait.area;
        for stem in branches
            #stem.xylem.state.kmax = stem.xylem.trait.area / trunk.xylem.trait.area * ks[3] * stem.xylem.trait.l / stem.xylem.trait.area;
            stem.xylem.trait.k_max = ks[3] * stem.xylem.trait.l / trunk.xylem.trait.area;
        end;
        for leaf in leaves
            leaf.xylem.trait.k_max = ks[4] / (can_str.trait.lai * sbulk.trait.area);
        end;
    end;

    # update Vcmax and Jmax TD if leaf shedding flag is not true
    if !spac.plant._leaf_shedded
        for leaf in leaves
            prescribe_ps_td!(config, leaf.photosystem.trait; t_clm = t_clm);
        end;
    end;

    #
    # canopy structure
    #
    # update LAI and leaf area if leaf shedding flag is not true
    if !isnothing(lai) && !spac.plant._leaf_shedded
        epslai = max(0, lai);
        can_str.trait.lai = epslai;
        can_str.trait.δlai = epslai .* ones(FT, n_layer) ./ n_layer;
        for irt in 1:n_layer
            ilf = n_layer - irt + 1;
            leaves[ilf].xylem.trait.area = sbulk.trait.area * can_str.trait.δlai[irt];
        end;
    end;

    # update CI
    if !isnothing(ci)
        can_str.trait.ci = ci;
    end;

    # update SAI
    if !isnothing(sai)
        can_str.trait.sai = sai;
        can_str.trait.δsai = sai .* ones(FT, n_layer) ./ n_layer;
    end;

    #
    # leaf temperature
    #
    # prescribe leaf temperature if leaf shedding flag is not true
    if !isnothing(t_leaf) && !spac.plant._leaf_shedded
        for leaf in leaves
            leaf.energy.s_aux.t = t_leaf;
        end;
    end;

    # update vcmax25 at the top layer (last element of leaves array because leaves are ordered from bottom to top) if leaf shedding flag is not true
    if !spac.plant._leaf_shedded
        prescribe_ps_traits!(leaves[end]; b6f = b6f, jmax = jmax, rd = rd, vcmax = vcmax, vpmax = vpmax);
        prescribe_ps_traits!(spac; vertical_expo = vertical_expo);
    end;

    # re-initialize leaf energy if LAI or t_leaf is updated (if leaf shedding flag is not true)
    if (!isnothing(lai) || !isnothing(t_leaf)) && !spac.plant._leaf_shedded
        for leaf in leaves
            initialize_energy_states!(leaf);
        end;
    end;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2024-Jul-23: add method to prescribe Vcmax, Jmax, b6f, and Rd for C3 models (C4 models pending)
#
#######################################################################################################################################################################################################
"""

    prescribe_ps_traits!(
                leaf::Leaf;
                b6f::Union{Nothing, Number} = nothing,
                jmax::Union{Nothing, Number} = nothing,
                rd::Union{Nothing, Number} = nothing,
                vcmax::Union{Nothing, Number} = nothing,
                vpmax::Union{Nothing, Number} = nothing)
    prescribe_ps_traits!(spac::BulkSPAC; vertical_expo::Union{Nothing, Number} = nothing)

Prescribe the photosynthetic traits for a single leaf, given
- `leaf` Leaf object
- `b6f` b6f content. Optional, default is nothing
- `jmax` Jmax25. Optional, default is nothing
- `rd` Dark respiration rate. Optional, default is nothing
- `vcmax` Vcmax25. Optional, default is nothing
- `vpmax` Vpmax25. Optional, default is nothing
- `spac` BulkSPAC object
- `vertical_expo` Exponential tuning factor to adjust Vcmax25. Optional, default is nothing

"""
function prescribe_ps_traits! end;

# Prescribe the variables for single leaf, suggest to do this only for top canopy
prescribe_ps_traits!(
            leaf::Leaf;
            b6f::Union{Nothing, Number} = nothing,
            jmax::Union{Nothing, Number} = nothing,
            rd::Union{Nothing, Number} = nothing,
            vcmax::Union{Nothing, Number} = nothing,
            vpmax::Union{Nothing, Number} = nothing) = prescribe_ps_traits!(leaf.photosystem.trait; b6f=b6f, jmax=jmax, rd=rd, vcmax=vcmax, vpmax=vpmax);

prescribe_ps_traits!(
            pst::Union{C3CLMTrait, C3FvCBTrait, C3VJPTrait};
            b6f::Union{Nothing, Number} = nothing,
            jmax::Union{Nothing, Number} = nothing,
            rd::Union{Nothing, Number} = nothing,
            vcmax::Union{Nothing, Number} = nothing,
            vpmax::Union{Nothing, Number} = nothing) = (
    if !isnothing(vcmax)
        pst.v_cmax25 = vcmax;

        # if jmax is not nothing
        if !isnothing(jmax)
            pst.j_max25 = jmax;
        else
            pst.j_max25 = vcmax * 1.6;
        end;

        # if rd is not nothing
        if !isnothing(rd)
            pst.r_d25 = rd;
        else
            pst.r_d25 = vcmax * 0.015;
        end;
    else
        # if jmax is not nothing
        if !isnothing(jmax)
            pst.j_max25 = jmax;
        end;

        # if rd is not nothing
        if !isnothing(rd)
            pst.r_d25 = rd;
        end;
    end;

    return nothing
);

prescribe_ps_traits!(
            pst::Union{C3CytoMinEtaTrait, C3CytoTrait, C3JBTrait};
            b6f::Union{Nothing, Number} = nothing,
            jmax::Union{Nothing, Number} = nothing,
            rd::Union{Nothing, Number} = nothing,
            vcmax::Union{Nothing, Number} = nothing,
            vpmax::Union{Nothing, Number} = nothing) = (
    if !isnothing(vcmax)
        pst.v_cmax25 = vcmax;

        # if b6f is not nothing
        if !isnothing(b6f)
            pst.b₆f = b6f;
        else
            pst.b₆f = vcmax * 0.0062;
        end;

        # if rd is not nothing
        if !isnothing(rd)
            pst.r_d25 = rd;
        else
            pst.r_d25 = vcmax * 0.015;
        end;
    else
        # if b6f is not nothing
        if !isnothing(b6f)
            pst.b₆f = b6f;
        end;

        # if rd is not nothing
        if !isnothing(rd)
            pst.r_d25 = rd;
        end;
    end;

    return nothing
);

prescribe_ps_traits!(
            pst::C4CLMTrait;
            b6f::Union{Nothing, Number} = nothing,
            jmax::Union{Nothing, Number} = nothing,
            rd::Union{Nothing, Number} = nothing,
            vcmax::Union{Nothing, Number} = nothing,
            vpmax::Union{Nothing, Number} = nothing) = (
    if !isnothing(vcmax)
        pst.v_cmax25 = vcmax;

        # if rd is not nothing
        if !isnothing(rd)
            pst.r_d25 = rd;
        else
            pst.r_d25 = vcmax * 0.015;
        end;
    else
        # if rd is not nothing
        if !isnothing(rd)
            pst.r_d25 = rd;
        end;
    end;

    return nothing
);

prescribe_ps_traits!(
            pst::C4VJPTrait;
            b6f::Union{Nothing, Number} = nothing,
            jmax::Union{Nothing, Number} = nothing,
            rd::Union{Nothing, Number} = nothing,
            vcmax::Union{Nothing, Number} = nothing,
            vpmax::Union{Nothing, Number} = nothing) = (
    if !isnothing(vcmax)
        pst.v_cmax25 = vcmax;
    end;

    if !isnothing(vpmax)
        pst.v_pmax25 = vpmax;

        # if rd is not nothing
        if !isnothing(rd)
            pst.r_d25 = rd;
        else
            pst.r_d25 = vpmax * 0.015;
        end;
    else
        # if rd is not nothing
        if !isnothing(rd)
            pst.r_d25 = rd;
        end;
    end;

    return nothing
);

# Method to apply the exponential tuning factor to Vcmax25...
prescribe_ps_traits!(spac::BulkSPAC; vertical_expo::Union{Nothing, Number} = nothing) = prescribe_ps_traits!(spac, spac.plant.leaves[end].photosystem.trait; vertical_expo=vertical_expo);

prescribe_ps_traits!(spac::BulkSPAC, pst::Union{C3CLMTrait, C3FvCBTrait, C3VJPTrait}; vertical_expo::Union{Nothing, Number} = nothing) = (
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # update vertical profiles
    for irt in 1:n_layer
        ilf = n_layer - irt + 1;
        ratio = isnothing(vertical_expo) ? 1 : exp(-vertical_expo * sum(can_str.trait.δlai[1:irt-1]));
        leaf = leaves[ilf];
        leaf.photosystem.trait.v_cmax25 = leaves[end].photosystem.trait.v_cmax25 * ratio;
        leaf.photosystem.trait.j_max25 = leaves[end].photosystem.trait.j_max25 * ratio;
        leaf.photosystem.trait.r_d25 = leaves[end].photosystem.trait.r_d25 * ratio;
    end;

    return nothing
);

prescribe_ps_traits!(spac::BulkSPAC, pst::Union{C3CytoMinEtaTrait, C3CytoTrait, C3JBTrait}; vertical_expo::Union{Nothing, Number} = nothing) = (
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # update vertical profiles
    for irt in 1:n_layer
        ilf = n_layer - irt + 1;
        ratio = isnothing(vertical_expo) ? 1 : exp(-vertical_expo * sum(can_str.trait.δlai[1:irt-1]));
        leaf = leaves[ilf];
        leaf.photosystem.trait.v_cmax25 = leaves[end].photosystem.trait.v_cmax25 * ratio;
        leaf.photosystem.trait.b₆f = leaves[end].photosystem.trait.b₆f * ratio;
        leaf.photosystem.trait.r_d25 = leaves[end].photosystem.trait.r_d25 * ratio;
    end;

    return nothing
);

prescribe_ps_traits!(spac::BulkSPAC, pst::C4CLMTrait; vertical_expo::Union{Nothing, Number} = nothing) = (
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # update vertical profiles
    for irt in 1:n_layer
        ilf = n_layer - irt + 1;
        ratio = isnothing(vertical_expo) ? 1 : exp(-vertical_expo * sum(can_str.trait.δlai[1:irt-1]));
        leaf = leaves[ilf];
        leaf.photosystem.trait.v_cmax25 = leaves[end].photosystem.trait.v_cmax25 * ratio;
        leaf.photosystem.trait.r_d25 = leaves[end].photosystem.trait.r_d25 * ratio;
    end;

    return nothing
);

prescribe_ps_traits!(spac::BulkSPAC, pst::C4VJPTrait; vertical_expo::Union{Nothing, Number} = nothing) = (
    can_str = spac.canopy.structure;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # update vertical profiles
    for irt in 1:n_layer
        ilf = n_layer - irt + 1;
        ratio = isnothing(vertical_expo) ? 1 : exp(-vertical_expo * sum(can_str.trait.δlai[1:irt-1]));
        leaf = leaves[ilf];
        leaf.photosystem.trait.v_cmax25 = leaves[end].photosystem.trait.v_cmax25 * ratio;
        leaf.photosystem.trait.v_pmax25 = leaves[end].photosystem.trait.v_pmax25 * ratio;
        leaf.photosystem.trait.r_d25 = leaves[end].photosystem.trait.r_d25 * ratio;
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2024-Jul-23: add method to prescribe Vcmax25 and Jmax25 TD temperature dependent
#
#######################################################################################################################################################################################################
"""

    prescribe_ps_td!(
                config::SPACConfiguration{FT},
                pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}, C3CytoMinEtaTrait{FT}, C3CytoTrait{FT}, C3JBTrait{FT}, C4CLMTrait{FT}, C4VJPTrait{FT}};
                t_clm::Union{Nothing, Number}) where {FT}

Prescribe the photosystem temperature dependence, given
- `config` Configuration for SPAC
- `pst` Photosystem trait type
- `t_clm` Moving average temperature to update Vcmax and Jmax temperature dependencies

"""
function prescribe_ps_td! end;

prescribe_ps_td!(config::SPACConfiguration{FT}, pst::Union{C3CLMTrait{FT}, C3FvCBTrait{FT}, C3VJPTrait{FT}}; t_clm::Union{Nothing, Number}) where {FT} = (
    (; T_CLM) = config;

    if !isnothing(t_clm)
        if T_CLM
            pst.TD_VCMAX.ΔSV = 668.39 - 1.07 * (t_clm - T₀(FT));
            pst.TD_JMAX.ΔSV = 659.70 - 0.75 * (t_clm - T₀(FT));
        end;
    end;

    return nothing
);

prescribe_ps_td!(config::SPACConfiguration{FT}, pst::Union{C3CytoMinEtaTrait{FT}, C3CytoTrait{FT}, C3JBTrait{FT}, C4CLMTrait{FT}, C4VJPTrait{FT}}; t_clm::Union{Nothing, Number}) where {FT} = (
    (; T_CLM) = config;

    if !isnothing(t_clm)
        if T_CLM
            pst.TD_VCMAX.ΔSV = 668.39 - 1.07 * (t_clm - T₀(FT));
        end;
    end;

    return nothing
);
