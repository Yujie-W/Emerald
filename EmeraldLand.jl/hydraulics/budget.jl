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
        ROOTS[_i].NS.energy.state.energy += ROOTS[_i].NS.energy.auxil.∂e∂t * δt;
        ROOTS[_i].NS.energy.auxil.cp = sum(ROOTS[_i].NS.xylem.state.v_storage) * CP_L_MOL(FT) + ROOTS[_i].NS.xylem.state.cp * ROOTS[_i].NS.xylem.state.area * ROOTS[_i].NS.xylem.state.l;
        ROOTS[_i].NS.energy.auxil.t = ROOTS[_i].NS.energy.state.energy / ROOTS[_i].NS.energy.auxil.cp;

        if isnan(ROOTS[_i].NS.energy.auxil.t) || isnan(ROOTS[_i].NS.energy.state.energy)
            @info "Debugging" ROOTS[_i].NS.energy.auxil.t ROOTS[_i].NS.energy.state.energy sum(ROOTS[_i].NS.xylem.state.v_storage);
            error("NaN detected when updating temperature for root");
        end;
    end;

    # update the temperature for trunk
    TRUNK.NS.energy.state.energy += TRUNK.NS.energy.auxil.∂e∂t * δt;
    TRUNK.NS.energy.auxil.cp = sum(TRUNK.NS.xylem.state.v_storage) * CP_L_MOL(FT) + TRUNK.NS.xylem.state.cp * TRUNK.NS.xylem.state.area * TRUNK.NS.xylem.state.l;
    TRUNK.NS.energy.auxil.t = TRUNK.NS.energy.state.energy / TRUNK.NS.energy.auxil.cp;

    if isnan(TRUNK.NS.energy.auxil.t) || isnan(TRUNK.NS.energy.state.energy)
        @info "Debugging" TRUNK.NS.energy.auxil.t TRUNK.NS.energy.state.energy sum(TRUNK.NS.xylem.state.v_storage);
        error("NaN detected when updating temperature for trunk");
    end;

    # update the temperature for branches and leaves
    for _i in 1:DIM_LAYER
        BRANCHES[_i].NS.energy.state.energy += BRANCHES[_i].NS.energy.auxil.∂e∂t * δt;
        BRANCHES[_i].NS.energy.auxil.cp = sum(BRANCHES[_i].NS.xylem.state.v_storage) * CP_L_MOL(FT) + BRANCHES[_i].NS.xylem.state.cp * BRANCHES[_i].NS.xylem.state.area * BRANCHES[_i].NS.xylem.state.l;
        BRANCHES[_i].NS.energy.auxil.t = BRANCHES[_i].NS.energy.state.energy / BRANCHES[_i].NS.energy.auxil.cp;
        LEAVES[_i].NS.energy.state.energy += LEAVES[_i].NS.energy.auxil.∂e∂t * δt;
        LEAVES[_i].NS.energy.auxil.cp = LEAVES[_i].NS.capacitor.state.v_storage * CP_L_MOL(FT) + LEAVES[_i].NS.xylem.state.cp * LEAVES[_i].NS.xylem.state.area * LEAVES[_i].NS.bio.state.lma * 10;
        LEAVES[_i].NS.energy.auxil.t = LEAVES[_i].NS.energy.state.energy / LEAVES[_i].NS.energy.auxil.cp;

        if isnan(BRANCHES[_i].NS.energy.auxil.t) || isnan(BRANCHES[_i].NS.energy.state.energy) || isnan(LEAVES[_i].NS.energy.auxil.t) || isnan(LEAVES[_i].NS.energy.state.energy)
            @info "Debugging" BRANCHES[_i].NS.energy.auxil.t BRANCHES[_i].NS.energy.state.energy;
            @info "Debugging" LEAVES[_i].NS.energy.auxil.t LEAVES[_i].NS.energy.state.energy;
            error("NaN detected when updating temperature for branch and leaf");
        end;
    end;

    return nothing
);
