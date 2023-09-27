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


"""

    plant_energy!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Compute the marginal energy increase in spac, given
- `config` SPAC configurations
- `spac` `MultiLayerSPAC` type SPAC

"""
plant_energy!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; DIM_LAYER, DIM_ROOT) = config;
    (; AIR, BRANCHES, CANOPY, LEAVES, LEAVES_INDEX, ROOTS, ROOTS_INDEX, SOIL, TRUNK) = spac;

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

        if isnan(ROOTS[_i].∂e∂t) || isnan(TRUNK.∂e∂t)
            @info "Debugging" ROOTS[_i].∂e∂t TRUNK.∂e∂t flow_in(ROOTS[_i]) flow_out(ROOTS[_i]) ROOTS[_i].t TRUNK.t SOIL.LAYERS[ROOTS_INDEX[_i]].t;
            error("NaN detected when computing energy budget for root and trunk");
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

        if isnan(BRANCHES[_i].∂e∂t) || isnan(LEAVES[_i].∂e∂t) || isnan(TRUNK.∂e∂t)
            @info "Debugging" BRANCHES[_i].∂e∂t LEAVES[_i].∂e∂t TRUNK.∂e∂t flow_in(BRANCHES[_i]) flow_out(BRANCHES[_i]) TRUNK.t BRANCHES[_i].t LEAVES[_i].t;
            error("NaN detected when computing energy budget for trunk, branch, and leaf");
        end;
    end;

    # loop through the leaves
    if CANOPY.lai == 0
        for _i in 1:DIM_LAYER
            LEAVES[_i].∂e∂t = 0;
        end;
    else
        for _i in 1:DIM_LAYER
            _g_be = FT(1.4) * FT(0.135) * sqrt(AIR[LEAVES_INDEX[_i]].wind / (FT(0.72) * LEAVES[_i].WIDTH));

            LEAVES[_i].∂e∂t  = 0;
            LEAVES[_i].∂e∂t += (CANOPY.RADIATION.r_net_sw[DIM_LAYER+1-_i] + CANOPY.RADIATION.r_net_lw[DIM_LAYER+1-_i]) / CANOPY.δlai[_i];
            LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * M_H₂O(FT) * latent_heat_vapor(LEAVES[_i].t);
            LEAVES[_i].∂e∂t -= flow_out(LEAVES[_i]) * CP_V_MOL(FT) * LEAVES[_i].t; # note here that CP_L_MOL is included in the latent_heat_vapor TD function
            LEAVES[_i].∂e∂t -= 2 * _g_be * CP_D_MOL(FT) * (LEAVES[_i].t - AIR[LEAVES_INDEX[_i]].t);

            if isnan(LEAVES[_i].∂e∂t)
                @info "Debugging" LEAVES[_i].∂e∂t CANOPY.RADIATION.r_net_sw[DIM_LAYER+1-_i] CANOPY.RADIATION.r_net_lw[DIM_LAYER+1-_i] CANOPY.δlai[_i];
                @info "Debugging" flow_out(LEAVES[_i]) latent_heat_vapor(LEAVES[_i].t) LEAVES[_i].t _g_be AIR[LEAVES_INDEX[_i]].t;
                error("NaN detected when computing energy budget for leaf");
            end;
        end;
    end;

    return nothing
);


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
        ROOTS[_i].e += ROOTS[_i].∂e∂t * δt;
        ROOTS[_i].t  = ROOTS[_i].e / (CP_L_MOL(FT) * sum(ROOTS[_i].HS.v_storage));

        if isnan(ROOTS[_i].t) || isnan(ROOTS[_i].e)
            @info "Debugging" ROOTS[_i].t ROOTS[_i].e sum(ROOTS[_i].HS.v_storage);
            error("NaN detected when updating temperature for root");
        end;
    end;

    # update the temperature for trunk
    TRUNK.e += TRUNK.∂e∂t * δt;
    TRUNK.t  = TRUNK.e / (CP_L_MOL(FT) * sum(TRUNK.HS.v_storage));

    if isnan(TRUNK.t) || isnan(TRUNK.e)
        @info "Debugging" TRUNK.t TRUNK.e sum(TRUNK.HS.v_storage);
        error("NaN detected when updating temperature for trunk");
    end;

    # update the temperature for branches and leaves
    for _i in 1:DIM_LAYER
        BRANCHES[_i].e += BRANCHES[_i].∂e∂t * δt;
        BRANCHES[_i].t  = BRANCHES[_i].e / (CP_L_MOL(FT) * sum(BRANCHES[_i].HS.v_storage));
        LEAVES[_i].e   += LEAVES[_i].∂e∂t * δt;
        LEAVES[_i].t    = LEAVES[_i].e / (LEAVES[_i].CP * LEAVES[_i].BIO.state.lma * 10 + CP_L_MOL(FT) * LEAVES[_i].HS.v_storage);

        if isnan(BRANCHES[_i].t) || isnan(BRANCHES[_i].e) || isnan(LEAVES[_i].t) || isnan(LEAVES[_i].e)
            @info "Debugging" BRANCHES[_i].t BRANCHES[_i].e sum(BRANCHES[_i].HS.v_storage);
            @info "Debugging" LEAVES[_i].t LEAVES[_i].e sum(LEAVES[_i].HS.v_storage);
            error("NaN detected when updating temperature for branch and leaf");
        end;
    end;

    return nothing
);
