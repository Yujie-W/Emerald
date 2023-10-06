#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-27: add function to initialize SPAC
#     2022-Jun-27: add leaf area controller to make sure soil and leaf areas are consistent with leaf area index
#     2023-Mar-27: initialize soil and leaf e as well (because T, SWC may be changed)
#     2023-Apr-13: add config to function call
#     2023-May-19: use δlai per canopy layer
#     2023-Jun-12: initialize soil trace gas as well
#     2023-Jun-13: update N₂ and O₂ based on soil water content
#     2023-Jun-13: add soil gas energy into soil e
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Sep-07: add ALLOW_SOIL_EVAPORATION check
#     2023-Sep-11: rename ALLOW_SOIL_EVAPORATION to ENABLE_SOIL_EVAPORATION
#
#######################################################################################################################################################################################################
"""

    initialize!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Initialize the SPAC, given
- `config` Configurations of spac model
- `spac` `MultiLayerSPAC` SPAC

"""
function initialize! end

initialize!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT} = (
    (; ENABLE_SOIL_EVAPORATION) = config;
    (; CANOPY, LEAVES, SOIL_BULK, SOILS) = spac;

    # make sure soil energy is correctly scaled with temperature and soil water content
    for soil in SOILS
        if ENABLE_SOIL_EVAPORATION
            δθ = max(0, soil.state.vc.Θ_SAT - soil.state.θ);
            rt = GAS_R(FT) * soil.auxil.t;
            soil.state.ns[3] = saturation_vapor_pressure(soil.auxil.t, soil.auxil.ψ * 1000000) * soil.auxil.δz * δθ / rt;
            soil.state.ns[4] = spac.AIR[1].P_AIR * 0.79 * soil.auxil.δz * δθ / rt;
            soil.state.ns[5] = spac.AIR[1].P_AIR * 0.209 * soil.auxil.δz * δθ / rt;
        end;
        cp_gas = (soil.state.ns[3] * CP_V_MOL(FT) + (soil.state.ns[1] + soil.state.ns[2] + soil.state.ns[4] + soil.state.ns[5]) * CP_D_MOL(FT)) / soil.auxil.δz;
        soil.auxil.cp = heat_capacitance(soil);
        soil.state.Σe = soil.auxil.cp * soil.auxil.t;
    end;

    # make sure leaf area index setup and energy are correct
    for i in eachindex(LEAVES)
        leaf = LEAVES[i];
        leaf.xylem.state.area = SOIL_BULK.state.area * CANOPY.δlai[i];
        initialize_energy_storage!(leaf);
    end;

    # initialize leaf level spectra
    plant_leaf_spectra!(config, spac);

    # initialize stomatal conductance
    stomatal_conductance!(spac, FT(0));

    return nothing
);
