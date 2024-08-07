# This file contains function to calculate energy budgets of the leaf

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2023-Sep-30: add function to calculate the energy flow of the leaf
#     2023-Nov-03: save latent and sensible heat fluxes in leaf energy auxil
#     2024-Feb-28: add LAI <= 0 control
#     2024-Jul-22: add option to enable chemical energy from photosynthesis or respiration
#
#######################################################################################################################################################################################################
"""

    leaf_energy_flows!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Calculate the energy flows of the leaf, given
- `config` `SPACConfiguration` type configuration
- `spac` `BulkSPAC` type SPAC

"""
function leaf_energy_flows! end;

leaf_energy_flows!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = leaf_energy_flows!(config, spac, spac.plant.leaves[1]);

leaf_energy_flows!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, ::CanopyLayer{FT}) where {FT} = (
    if spac.canopy.structure.trait.lai <= 0
        return nothing
    end;

    (; ENABLE_CHEMICAL_ENERGY) = config;

    # run the energy budget for each leaf layer only if LAI > 0
    airs = spac.airs;
    branches = spac.plant.branches;
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    lindex = spac.plant.leaves_index;
    sbulk = spac.soil_bulk;
    n_layer = length(leaves);

    # the total energy change of the leaf is the sum of
    #     the energy of the flow from the stem
    #     the energy of the flow from the air
    #     the sensible heat from the air
    #     the net radiation energy from shortwave and longwave radiation

    for irt in 1:n_layer
        ilf = n_layer + 1 - irt;
        stem = branches[ilf];
        leaf = leaves[ilf];
        air = airs[lindex[ilf]];

        # if flow in is positive, then energy flow is positive for leaf
        f_i = flow_in(leaf);
        if f_i >= 0
            leaf.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * stem.energy.s_aux.t;
        else
            leaf.energy.auxil.∂e∂t += f_i * CP_L_MOL(FT) * leaf.energy.s_aux.t;
        end;

        # if flow out is positive, then energy flow is negative for leaf
        # if f_o is positive, then the leaf is losing water to the air
        # note here that CP_L_MOL is included in the latent_heat_vapor TD function
        # so here we only need to calculate the heat mass flow from water vapor in gas phase
        f_o = flow_out(leaf);
        leaf.energy.auxil.∂e∂t -= f_o * M_H₂O(FT) * latent_heat_vapor(leaf.energy.s_aux.t);
        le = f_o >= 0 ? f_o * CP_V_MOL(FT) * leaf.energy.s_aux.t : f_o * CP_V_MOL(FT) * air.s_aux.t;
        leaf.energy.auxil.∂e∂t -= le;
        leaf.energy.auxil.∂e∂t_le = le;

        # add the sensible heat flux from the leaf to air (to total leaf area)
        g_be = FT(1.35) * leaf.flux.auxil.g_CO₂_b;
        sh = 2 * g_be * CP_D_MOL(FT) * (leaf.energy.s_aux.t - air.s_aux.t) * leaf.xylem.trait.area;
        leaf.energy.auxil.∂e∂t -= sh;
        leaf.energy.auxil.∂e∂t_sh = sh;

        # add the net radiation energy to the leaf (to total leaf area)
        leaf.energy.auxil.∂e∂t += (canopy.sun_geometry.auxil.r_net_sw_leaf[irt] + canopy.structure.auxil.r_net_lw_leaf[irt]) * sbulk.trait.area;

        # remove the chemical energy from the leaf
        if ENABLE_CHEMICAL_ENERGY
            an_layer = leaves[ilf].flux.auxil.a_n_mean * leaf.xylem.trait.area;
            leaf.energy.auxil.∂e∂t -= an_layer * FT(1e-6) / 6 * GLUCOSE(FT);
        end;
    end;

    return nothing
);
