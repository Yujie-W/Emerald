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
                config::SPACConfig{FT},
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
                t_acclim::Union{Number,Nothing} = nothing,
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
- `t_acclim` Moving average temperature to update Vcmax and Jmax temperature dependencies. Optional, default is nothing
- `t_leaf` Leaf temperature. Optional, default is nothing
- `vcmax` Vcmax25 at the top of canopy. Optional, default is nothing
- `vertical_expo` Exponential tuning factor to adjust Vcmax25. Optional, default is nothing
- `vpmax` Vpmax25 at the top of canopy. Optional, default is nothing

"""
function prescribe_traits!(
            config::SPACConfig{FT},
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
            t_acclim::Union{Number,Nothing} = nothing,
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
        prescribe_ps_td!(config; t_acclim = t_acclim);
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
            leaf.energy.auxil.t = t_leaf;
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
