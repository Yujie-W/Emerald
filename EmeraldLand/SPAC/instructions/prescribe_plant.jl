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
#
#######################################################################################################################################################################################################
"""

    prescribe_traits!(
                config::SPACConfiguration{FT},
                spac::BulkSPAC{FT};
                cab::Union{Number,Nothing} = nothing,
                car::Union{Number,Nothing} = nothing,
                ci::Union{Number,Nothing} = nothing,
                kmax::Union{Number,Tuple,Nothing} = nothing,
                lai::Union{Number,Nothing} = nothing,
                sai::Union{Number,Nothing} = nothing,
                swcs::Union{Tuple,Nothing} = nothing,
                t_clm::Union{Number,Nothing} = nothing,
                t_leaf::Union{Number,Nothing} = nothing,
                t_soils::Union{Tuple,Nothing} = nothing,
                vcmax::Union{Number,Nothing} = nothing,
                vcmax_expo::Union{Number,Nothing} = nothing) where {FT}

Update the physiological parameters of the SPAC, given
- `spac` Soil plant air continuum
- `config` Configuration for `BulkSPAC`
- `cab` Chlorophyll content. Optional, default is nothing
- `car` Carotenoid content. Optional, default is nothing
- `ci` Clumping index. Optional, default is nothing
- `kmax` Maximum hydraulic conductance. Optional, default is nothing
- `lai` Leaf area index. Optional, default is nothing
- `sai` Stem area index. Optional, default is nothing
- `t_clm` Moving average temperature to update Vcmax and Jmax temperature dependencies. Optional, default is nothing
- `t_leaf` Leaf temperature. Optional, default is nothing
- `vcmax` Vcmax25 at the top of canopy. Optional, default is nothing
- `vcmax_expo` Exponential tuning factor to adjust Vcmax25. Optional, default is nothing

"""
function prescribe_traits!(
            config::SPACConfiguration{FT},
            spac::BulkSPAC{FT};
            cab::Union{Number,Nothing} = nothing,
            car::Union{Number,Nothing} = nothing,
            ci::Union{Number,Nothing} = nothing,
            kmax::Union{Number,Tuple,Nothing} = nothing,
            lai::Union{Number,Nothing} = nothing,
            sai::Union{Number,Nothing} = nothing,
            t_clm::Union{Number,Nothing} = nothing,
            t_leaf::Union{Number,Nothing} = nothing,
            vcmax::Union{Number,Nothing} = nothing,
            vcmax_expo::Union{Number,Nothing} = nothing
) where {FT}
    (; T_CLM) = config;
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
    # update chlorophyll and carotenoid contents
    if !isnothing(cab)
        for leaf in leaves
            leaf.bio.trait.cab = cab;
        end;
    end;

    if !isnothing(car)
        for leaf in leaves
            leaf.bio.trait.car = car;
        end;
    end;

    # update kmax
    if !isnothing(kmax)
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

    # update Vcmax and Jmax TD
    if !isnothing(t_clm)
        for leaf in leaves
            if T_CLM
                leaf.photosystem.trait.TD_VCMAX.ΔSV = 668.39 - 1.07 * (t_clm - T₀(FT));
                leaf.photosystem.trait.TD_JMAX.ΔSV = 659.70 - 0.75 * (t_clm - T₀(FT));
            end;
        end;
    end;

    # update vcmax25 at the top layer (last element of leaves array because leaves are ordered from bottom to top)
    if !isnothing(vcmax)
        leaves[end].photosystem.trait.v_cmax25 = vcmax;
    end;

    #
    # canopy structure
    #
    # update LAI and leaf area
    if !isnothing(lai)
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
    # prescribe leaf temperature
    if !isnothing(t_leaf)
        for leaf in leaves
            leaf.energy.s_aux.t = t_leaf;
        end;
    end;

    #
    # parameters related to changes in any of multiple traits
    #
    # update vertical vcmax profile
    if !isnothing(vcmax) || !isnothing(lai)
        for irt in 1:n_layer
            ilf = n_layer - irt + 1;
            ratio = isnothing(vcmax_expo) ? 1 : exp(-vcmax_expo * sum(can_str.trait.δlai[1:irt-1]));
            leaf = leaves[ilf];
            if typeof(leaf.photosystem.trait) isa C3CytoTrait
                leaf.photosystem.trait.v_cmax25 = leaves[end].photosystem.trait.v_cmax25 * ratio;
                leaf.photosystem.trait.b₆f = leaves[end].photosystem.trait.v_cmax25 * 7 / 300 * ratio;
                leaf.photosystem.trait.r_d25 = leaves[end].photosystem.trait.v_cmax25 * 0.015 * ratio;
            elseif leaf.photosystem.trait isa C3VJPTrait
                leaf.photosystem.trait.v_cmax25 = leaves[end].photosystem.trait.v_cmax25 * ratio;
                leaf.photosystem.trait.j_max25 = leaves[end].photosystem.trait.v_cmax25 * 1.67 * ratio;
                leaf.photosystem.trait.r_d25 = leaves[end].photosystem.trait.v_cmax25 * 0.015 * ratio;
            else
                error("Vcmax profile is only available for C3CytoTrait and C3VJPTrait.");
            end;
        end;
    end;

    # re-initialize leaf energy if LAI or t_leaf is updated
    if !isnothing(lai) || !isnothing(t_leaf)
        for leaf in leaves
            initialize_energy_states!(leaf);
        end;
    end;

    return nothing
end;
