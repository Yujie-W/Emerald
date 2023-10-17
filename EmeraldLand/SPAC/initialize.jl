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
#     2023-Oct-07: add 0.01 to the water vapor volume per soil layer
#     2023-Oct-09: add root and stem initialization in the initialization of SPAC
#     2023-Oct-17: update step and subtep auxils during initialization
#
#######################################################################################################################################################################################################
"""

    initialize!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Initialize the SPAC, given
- `config` Configurations of spac model
- `spac` `BulkSPAC` SPAC

"""
function initialize! end;

initialize!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    airs = spac.airs;
    soils = spac.soils;
    soil_bulk = spac.soil_bulk;
    canopy = spac.canopy;

    roots = spac.plant.roots;
    trunk = spac.plant.trunk;
    branches = spac.plant.branches;
    leaves = spac.plant.leaves;

    # make sure soil energy is correctly scaled with temperature and soil water content
    for soil in soils
        initialize_struct!(soil, airs[1]);
    end;

    # make sure the root energy is correctly scaled with temperature
    for root in roots
        initialize_struct!(root);
    end;

    # make sure the stem energy is correctly scaled with temperature
    initialize_struct!(trunk);
    for stem in branches
        initialize_struct!(stem);
    end;

    # make sure leaf area index setup and energy are correct
    for i in eachindex(leaves)
        leaves[i].xylem.state.area = soil_bulk.state.area * canopy.structure.state.δlai[i];
        initialize_struct!(leaves[i]);
    end;

    # make sure air layers are correctly initialized
    for air in airs
        initialize_struct!(air);
    end;

    # make sure tha auxilary variables are correctly initialized
    update_step_auxils!(spac);
    update_substep_auxils!(spac);

    # initialize leaf level spectra
    plant_leaf_spectra!(config, spac);

    # initialize stomatal conductance
    stomatal_conductance!(spac, FT(0));

    return nothing
);
