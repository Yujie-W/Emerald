#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jul-15: add method for MultiLayerSPAC
#     2022-Jul-15: rename function to plant_energy! to be more accurate (ready to add other organs other leaf)
#     2022-Jul-15: add root, trunk, branch energy budgets
#     2022-Jul-26: add leaf LMA to the denominator
#     2022-Jul-27: fix the NaN temperature issue related to grass
#     2022-Jul-27: fix the unit issue related to CP and mol to kg
#     2022-Jul-27: rescale leaf energy absorption from per groud area to per leaf area
#     2023-Mar-27: fix a typo in 1:DIM_ROOT (was typed as DIM_ROOT)
#     2023-Mar-27: redo the energy balance related to water transport by make a statement of whether the flow is positive or negative
#     2023-May-19: use δlai per canopy layer
#     2023-Aug-27: add nan check
#     2023-Sep-07: bug fix for the energy budget of leaf as CP_V and CP_L are accounted for in the latent heat of vaporization function
#
#######################################################################################################################################################################################################
"""
This function has two major functionalities:
- Compute marginal energy increase in each organ
- Update the temperature in each organ when time step provided

"""
function plant_energy! end


# for _i in 1:DIM_LAYER
#     _g_be = FT(1.4) * FT(0.135) * sqrt(AIR[LEAVES_INDEX[_i]].wind / (FT(0.72) * LEAVES[_i].WIDTH));
#     LEAVES[_i].NS.energy.auxil.∂e∂t += (CANOPY.RADIATION.r_net_sw[DIM_LAYER+1-_i] + CANOPY.RADIATION.r_net_lw[DIM_LAYER+1-_i]) / CANOPY.δlai[_i];
#     LEAVES[_i].NS.energy.auxil.∂e∂t -= 2 * _g_be * CP_D_MOL(FT) * (LEAVES[_i].NS.energy.auxil.t - AIR[LEAVES_INDEX[_i]].t);


"""

    plant_energy!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Compute the marginal energy increase in spac, given
- `config` SPAC configurations
- `spac` `MultiLayerSPAC` type SPAC
- `δt` Time step

"""
plant_energy!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}, δt::FT) where {FT} = (
    (; DIM_LAYER, DIM_ROOT) = config;
    (; BRANCHES, LEAVES, ROOTS, TRUNK) = spac;

    # update the temperature for roots
    for _i in 1:DIM_ROOT
        ROOTS[_i].energy.state.Σe += ROOTS[_i].energy.auxil.∂e∂t * δt;
    end;

    # update the temperature for trunk
    TRUNK.energy.state.Σe += TRUNK.energy.auxil.∂e∂t * δt;

    # update the temperature for branches and leaves
    for _i in 1:DIM_LAYER
        BRANCHES[_i].energy.state.Σe += BRANCHES[_i].energy.auxil.∂e∂t * δt;
        LEAVES[_i].NS.energy.state.Σe += LEAVES[_i].NS.energy.auxil.∂e∂t * δt;
    end;

    return nothing
);
