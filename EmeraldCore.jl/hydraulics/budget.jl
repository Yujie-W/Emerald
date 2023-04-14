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
#     2023-Mar-27: rodo the energy balance related to water transport by make a statement of whether the flow is positive or negative
#
#######################################################################################################################################################################################################
"""
This function has two major functionalities:
- Compute marginal energy increase in each organ
- Update the temperature in each organ when time step provided

"""
function plant_energy! end


"""

    plant_energy!(spac::MultiLayerSPAC{FT}) where {FT}

Compute the marginal energy increase in spac, given
- `spac` `MultiLayerSPAC` type SPAC

"""
plant_energy!(spac::MultiLayerSPAC{FT}) where {FT} = (
    (; AIR, BRANCHES, CANOPY, DIM_LAYER, DIM_ROOT, LEAVES, LEAVES_INDEX, ROOTS, ROOTS_INDEX, SOIL, TRUNK) = spac;

    # loop through the roots
    TRUNK.∂e∂t = 0;
    for _i in 1:DIM_ROOT
        ROOTS[_i].∂e∂t = 0;
        if flow_in(ROOTS[_i]) >= 0
            ROOTS[_i].∂e∂t += flow_in(ROOTS[_i]) * CP_L_MOL(FT) * SOIL.LAYERS[ROOTS_INDEX[_i]].t;
        else
            ROOTS[_i].∂e∂t += flow_in(ROOTS[_i]) * CP_L_MOL(FT) * ROOTS[_i].t;
        end;
        if flow_out(ROOTS[_i]) >= 0
            ROOTS[_i].∂e∂t -= flow_out(ROOTS[_i]) * CP_L_MOL(FT) * ROOTS[_i].t;
            TRUNK.∂e∂t     += flow_out(ROOTS[_i]) * CP_L_MOL(FT) * ROOTS[_i].t;
        else
            ROOTS[_i].∂e∂t -= flow_out(ROOTS[_i]) * CP_L_MOL(FT) * TRUNK.t;
            TRUNK.∂e∂t     += flow_out(ROOTS[_i]) * CP_L_MOL(FT) * TRUNK.t;
        end;
    end;

    # loop through the branches
    for _i in 1:DIM_LAYER
        BRANCHES[_i].∂e∂t = 0;
        LEAVES[_i].∂e∂t   = 0;
        if flow_in(BRANCHES[_i]) >= 0
            TRUNK.∂e∂t        -= flow_in(BRANCHES[_i]) * CP_L_MOL(FT) * TRUNK.t;
            BRANCHES[_i].∂e∂t += flow_in(BRANCHES[_i]) * CP_L_MOL(FT) * TRUNK.t;
        else
            TRUNK.∂e∂t        -= flow_in(BRANCHES[_i]) * CP_L_MOL(FT) * BRANCHES[_i].t;
            BRANCHES[_i].∂e∂t += flow_in(BRANCHES[_i]) * CP_L_MOL(FT) * BRANCHES[_i].t;
        end;
        if flow_out(BRANCHES[_i]) >= 0
            BRANCHES[_i].∂e∂t -= flow_out(BRANCHES[_i]) * CP_L_MOL(FT) * BRANCHES[_i].t;
            LEAVES[_i].∂e∂t   += flow_in(LEAVES[_i]) * CP_L_MOL(FT) * BRANCHES[_i].t;
        else
            BRANCHES[_i].∂e∂t -= flow_out(BRANCHES[_i]) * CP_L_MOL(FT) * LEAVES[_i].t;
            LEAVES[_i].∂e∂t   += flow_in(LEAVES[_i]) * CP_L_MOL(FT) * LEAVES[_i].t;
        end;
    end;

    # loop through the leaves
    if CANOPY.lai == 0
        for _i in 1:DIM_LAYER
            LEAVES[_i].∂e∂t  = 0;
        end;
    else
        for _i in 1:DIM_LAYER
            _g_be = FT(1.4) * FT(0.135) * sqrt(AIR[LEAVES_INDEX[_i]].wind / (FT(0.72) * LEAVES[_i].WIDTH));

            LEAVES[_i].∂e∂t  = 0;
            LEAVES[_i].∂e∂t += (CANOPY.RADIATION.r_net_sw[DIM_LAYER+1-_i] + CANOPY.RADIATION.r_net_lw[DIM_LAYER+1-_i]) / (CANOPY.lai / DIM_LAYER);
            LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * M_H₂O(FT) * latent_heat_vapor(LEAVES[_i].t);
            LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * CP_L_MOL(FT) * LEAVES[_i].t;
            LEAVES[_i].∂e∂t -= 2 * _g_be * CP_D_MOL(FT) * (LEAVES[_i].t - AIR[LEAVES_INDEX[_i]].t);
        end;
    end;

    return nothing
);


"""

    plant_energy!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT}

Compute the marginal energy increase in spac, given
- `spac` `MultiLayerSPAC` type SPAC
- `δt` Time step

"""
plant_energy!(spac::MultiLayerSPAC{FT}, δt::FT) where {FT} = (
    (; BRANCHES, DIM_LAYER, DIM_ROOT, LEAVES, ROOTS, TRUNK) = spac;

    # update the temperature for roots
    for _i in 1:DIM_ROOT
        ROOTS[_i].e += ROOTS[_i].∂e∂t * δt;
        ROOTS[_i].t  = ROOTS[_i].e / (CP_L_MOL(FT) * sum(ROOTS[_i].HS.v_storage));
    end;

    # update the temperature for trunk
    TRUNK.e += TRUNK.∂e∂t * δt;
    TRUNK.t  = TRUNK.e / (CP_L_MOL(FT) * sum(TRUNK.HS.v_storage));

    # update the temperature for branches and leaves
    for _i in 1:DIM_LAYER
        BRANCHES[_i].e += BRANCHES[_i].∂e∂t * δt;
        BRANCHES[_i].t  = BRANCHES[_i].e / (CP_L_MOL(FT) * sum(BRANCHES[_i].HS.v_storage));
        LEAVES[_i].e   += LEAVES[_i].∂e∂t * δt;
        LEAVES[_i].t    = LEAVES[_i].e / (LEAVES[_i].CP * LEAVES[_i].BIO.lma * 10 + CP_L_MOL(FT) * LEAVES[_i].HS.v_storage);
    end;

    return nothing
);
