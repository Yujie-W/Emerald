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
#     2024-Jul-24: add leaf shedded flag to LAI prescription
#     2024-Aug-06: add leaf regrow flag to LAI prescription
#     2024-Aug-29: use carbon pool to update LAI (when LAI increases)
#     2024-Sep-03: make sure to update leaf asap as well when LAI is updated
#     2024-Sep-04: when lai_diff > 0, make sure carbon pool is not immediately used up and recover leaf xylem hydraulic system
#     2024-Sep-07: improve ci prescription to account for angular dependency
#
#######################################################################################################################################################################################################
"""

    prescribe_traits!(
                config::SPACConfiguration{FT},
                spac::BulkSPAC{FT};
                b6f::Union{Number,Nothing} = nothing,
                cab::Union{Number,Nothing} = nothing,
                car::Union{Number,Nothing} = nothing,
                ci::Union{Number,Vector,Nothing} = nothing,
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
            ci::Union{Number,Vector,Nothing} = nothing,
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
            root.xylem.trait.k_max = ks[1] * root.xylem.trait.l / trunk.xylem.state.asap;
        end;
        trunk.xylem.trait.k_max = ks[2] * trunk.xylem.trait.l / trunk.xylem.state.asap;
        for stem in branches
            #stem.xylem.state.kmax = stem.xylem.trait.area / trunk.xylem.trait.area * ks[3] * stem.xylem.trait.l / stem.xylem.trait.area;
            stem.xylem.trait.k_max = ks[3] * stem.xylem.trait.l / trunk.xylem.state.asap;
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
    # update LAI and leaf area if leaf shedding flag is not true or regrow flag is true
    # clear the legacy of leaves if regrow flag is true
    # TODO: use shed_leaves! and grow_leaves! functions in the future
    if !isnothing(lai)
        # if lai is not 0, grow new leaves is allowed
        lai_0 = can_str.trait.lai;
        lai_diff = lai - lai_0;
        c_demand = lai_diff * sbulk.trait.area * spac.plant.leaves[1].bio.trait.lma * 10000 / 30;
        c_allocable = spac.plant.pool.c_pool - spac.plant.pool.c_pool_min;
        if !spac.plant._leaf_shedded
            if lai_diff > 0
                if c_allocable <= 0
                    lai_diff = 0;
                elseif c_demand > c_allocable
                    lai_diff *= c_allocable / c_demand;
                end;
            end;
        elseif spac.plant._leaf_regrow
            # lai_diff > 0 for sure
            if 0 < c_demand <= c_allocable
                nothing
            # c_demand > c_allocable
            elseif c_allocable <= spac.plant.pool.c_pool_min
                c_actual = spac.plant.pool.c_pool / 2;
                lai_diff *= c_actual / c_demand;
            else # c_allocable > spac.plant.pool.c_pool_min
                lai_diff *= c_allocable / c_demand;
            end;
        end;

        if !spac.plant._leaf_shedded || spac.plant._leaf_regrow
            # update the leaf area
            can_str.trait.lai = lai_0 + lai_diff;
            can_str.trait.δlai = can_str.trait.lai .* ones(FT, n_layer) ./ n_layer;
            for irt in 1:n_layer
                ilf = n_layer - irt + 1;
                leaves[ilf].xylem.trait.area = sbulk.trait.area * can_str.trait.δlai[irt];
                leaves[ilf].xylem.state.asap = leaves[ilf].xylem.trait.area;
            end;

            # if lai_diff is positive, remove the energy from the carbon pool
            if lai_diff > 0
                c_mol = lai_diff * sbulk.trait.area * spac.plant.leaves[1].bio.trait.lma * 10000 / 30;
                spac.plant.pool.c_pool -= c_mol;
            end;

            # reset the flags and clear the legacy of leaves
            if lai_diff > 0
                spac.plant._leaf_regrow = false;
                spac.plant._leaf_shedded = false;
                for l in leaves
                    xylem_recovery!(l.xylem, lai_0, lai_diff);
                end;
            end;
        end;
    end;

    # update CI
    if !isnothing(ci)
        if ci isa Number
            can_str.trait.ci.ci_0 = ci;
            can_str.trait.ci.ci_1 = 0;
        else
            can_str.trait.ci.ci_0 = ci[0];
            can_str.trait.ci.ci_1 = ci[1];
        end;
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
#     2024-Aug-01: use GeneralC3Trait and GeneralC4Trait
#     2024-Aug-09: typo fix for GeneralC4Trait
#
#######################################################################################################################################################################################################
"""

    prescribe_ps_traits!(
                leaf::Union{CanopyLayer, Leaf};
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
            leaf::Union{CanopyLayer, Leaf};
            b6f::Union{Nothing, Number} = nothing,
            jmax::Union{Nothing, Number} = nothing,
            rd::Union{Nothing, Number} = nothing,
            vcmax::Union{Nothing, Number} = nothing,
            vpmax::Union{Nothing, Number} = nothing) = prescribe_ps_traits!(leaf.photosystem.trait; b6f = b6f, jmax = jmax, rd = rd, vcmax = vcmax, vpmax = vpmax);

prescribe_ps_traits!(
            pst::Union{GeneralC3Trait, GeneralC4Trait};
            b6f::Union{Nothing, Number} = nothing,
            jmax::Union{Nothing, Number} = nothing,
            rd::Union{Nothing, Number} = nothing,
            vcmax::Union{Nothing, Number} = nothing,
            vpmax::Union{Nothing, Number} = nothing) = prescribe_ps_traits!(pst, pst.ACM, pst.AJM, pst.APM; b6f = b6f, jmax = jmax, rd = rd, vcmax = vcmax, vpmax = vpmax);

prescribe_ps_traits!(
            pst::GeneralC3Trait,
            acm::AcMethodC3VcmaxPi,
            ajm::AjMethodC3JmaxPi,
            apm::Union{ApMethodC3Inf, ApMethodC3Vcmax};
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
            pst.j_max25 = vcmax * 1.64;
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
            pst::GeneralC3Trait,
            acm::AcMethodC3VcmaxPi,
            ajm::AjMethodC3VqmaxPi,
            apm::Union{ApMethodC3Inf, ApMethodC3Vcmax};
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
            pst.b₆f = vcmax * 0.0066;
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
            pst::GeneralC4Trait,
            acm::AcMethodC4Vcmax,
            ajm::AjMethodC4JPSII,
            apm::ApMethodC4VcmaxPi;
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
            pst::GeneralC4Trait,
            acm::AcMethodC4Vcmax,
            ajm::AjMethodC4JPSII,
            apm::ApMethodC4VpmaxPi;
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
prescribe_ps_traits!(
            spac::BulkSPAC;
            vertical_expo::Union{Nothing, Number} = nothing) = prescribe_ps_traits!(spac, spac.plant.leaves[end].photosystem.trait; vertical_expo = vertical_expo);

prescribe_ps_traits!(
            spac::BulkSPAC,
            pst::Union{GeneralC3Trait, GeneralC4Trait};
            vertical_expo::Union{Nothing, Number} = nothing) = prescribe_ps_traits!(spac, pst.ACM, pst.AJM, pst.APM; vertical_expo = vertical_expo);

prescribe_ps_traits!(
            spac::BulkSPAC,
            acm::AcMethodC3VcmaxPi,
            ajm::AjMethodC3JmaxPi,
            apm::Union{ApMethodC3Inf, ApMethodC3Vcmax};
            vertical_expo::Union{Nothing, Number} = nothing) = (
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

prescribe_ps_traits!(
            spac::BulkSPAC,
            acm::AcMethodC3VcmaxPi,
            ajm::AjMethodC3VqmaxPi,
            apm::Union{ApMethodC3Inf, ApMethodC3Vcmax};
            vertical_expo::Union{Nothing, Number} = nothing) = (
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

prescribe_ps_traits!(
            spac::BulkSPAC,
            acm::AcMethodC4Vcmax,
            ajm::AjMethodC4JPSII,
            apm::ApMethodC4VcmaxPi;
            vertical_expo::Union{Nothing, Number} = nothing) = (
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

prescribe_ps_traits!(
            spac::BulkSPAC,
            acm::AcMethodC4Vcmax,
            ajm::AjMethodC4JPSII,
            apm::ApMethodC4VpmaxPi;
            vertical_expo::Union{Nothing, Number} = nothing) = (
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
#     2024-Aug-01: use GeneralC3Trait and GeneralC4Trait
#     2024-Aug-01: fix the dispatching issue with GeneralC4Trait models
#
#######################################################################################################################################################################################################
"""

    prescribe_ps_td!(
                config::SPACConfiguration{FT},
                pst::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}};
                t_clm::Union{Nothing, Number}) where {FT}

Prescribe the photosystem temperature dependence, given
- `config` Configuration for SPAC
- `pst` Photosystem trait type
- `t_clm` Moving average temperature to update Vcmax and Jmax temperature dependencies

"""
function prescribe_ps_td! end;

prescribe_ps_td!(
            config::SPACConfiguration{FT},
            pst::GeneralC3Trait{FT};
            t_clm::Union{Nothing, Number}) where {FT} = prescribe_ps_td!(config, pst.TD_VCMAX, pst.TD_JMAX; t_clm = t_clm);

prescribe_ps_td!(
            config::SPACConfiguration{FT},
            vtd::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}},
            jtd::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}};
            t_clm::Union{Nothing, Number}) where {FT} = (
    (; T_CLM) = config;

    # TODO: double check if there is memory allocation in this method
    if !isnothing(t_clm) && T_CLM
        if vtd isa ArrheniusPeak || vtd isa Q10Peak
            vtd.ΔSV = 668.39 - 1.07 * (t_clm - T₀(FT));
        end;
        if jtd isa ArrheniusPeak || jtd isa Q10Peak
            jtd.ΔSV = 659.70 - 0.75 * (t_clm - T₀(FT));
        end;
    end;

    return nothing
);

prescribe_ps_td!(
            config::SPACConfiguration{FT},
            pst::GeneralC4Trait{FT};
            t_clm::Union{Nothing, Number}) where {FT} = prescribe_ps_td!(config, pst.TD_VCMAX; t_clm = t_clm);

prescribe_ps_td!(
            config::SPACConfiguration{FT},
            vtd::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}, Q10Peak{FT}, Q10PeakHT{FT}, Q10PeakLTHT{FT}};
            t_clm::Union{Nothing, Number}) where {FT} = (
    (; T_CLM) = config;

    if !isnothing(t_clm) && T_CLM
        if vtd isa ArrheniusPeak || vtd isa Q10Peak
            vtd.ΔSV = 668.39 - 1.07 * (t_clm - T₀(FT));
        end;
    end;

    return nothing
);
