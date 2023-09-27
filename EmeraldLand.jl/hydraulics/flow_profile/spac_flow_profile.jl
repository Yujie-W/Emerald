#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-27: migrate function to new version
#     2022-May-27: rename the functions (flow_profile! and update_PVF!) to xylem_flow_profile!
#     2022-May-27: add method to set up root flow rate at steady state mode (this function is used to solve for steady state solution)
#     2022-May-27: add method to set up root flow rate at non-steady state mode (this function is used to solve for steady state solution)
#     2022-May-27: add method for leaf, root, and stem at steady state mode
#     2022-May-27: add method for leaf at non-steady state mode
#     2022-May-27: add method for root and stem at non-steady state mode
#     2022-May-27: add method for leaf, root, and stem hydraulic system at steady and non-steady state mode (for dispatching purpose)
#     2022-May-27: add method for leaf, root, and stem organ at steady and non-steady state mode (for dispatching purpose)
#     2022-May-27: add method to solve root flow rate partition at both steady and non-steady state modes
#     2022-May-27: add method for MonoElementSPAC (blank)
#     2022-May-31: remove hydraulic system from input variables, thus supporting leaf and stem
#     2022-May-31: use reformulate methods for setting up flow rate
#     2022-May-31: set up the flow rate profile using the network
#     2022-May-31: add method for MonoGrassSPAC
#     2022-May-31: add method for MonoPalmSPAC
#     2022-May-31: add method for MonoTreeSPAC
#     2022-Jun-29: rename SPAC to ML*SPAC to be more accurate
#     2022-Jun-30: add support to Leaves2D
#     2022-Jun-30: add method for Leaves1D
#     2022-Jul-12: add method to update leaf hydraulic flow rates per canopy layer based on stomatal conductance
#     2022-Oct-20: use add SoilLayer to function variables, because of the removal of SH from RootHydraulics
#     2022-Oct-20: fix a bug in flow profile counter (does not impact simulation)
#     2022-Oct-21: add a second solver to fix the case when root_pk does not work
#     2023-Mar-28: if root is disconnected, do not update its flow profile
#     2023-Mar-28: if no root is connected to soil, set all flow to 0
#     2023-Jun-16: compute saturated vapor pressure based on water water potential
#     2023-Aug-23: add configuration to enable/disable leaf condensation
#     2023-Sep-11: add config to the variable list
#     2023-Sep-11: rename methods to different functions to be more logical
#     2023-Sep-12: rename function to spac_flow_profile!
#
#######################################################################################################################################################################################################
"""

    spac_flow_profile!(config::SPACConfiguration{FT}, spac::MonoElementSPAC{FT}, Δt::FT) where {FT}

Update flow profiles for the soil-plant-air continuum (set up leaf flow rate from stomatal conductance first), given
- `config` `SPACConfiguration` type struct
- `spac` `MonoElementSPAC` or `MultiLayerSPAC` type SPAC system
- `Δt` Time step length

"""
function spac_flow_profile! end

spac_flow_profile!(config::SPACConfiguration{FT}, spac::Union{MonoElementSPAC{FT},MultiLayerSPAC{FT}}, Δt::FT) where {FT} = (
    # 1. update the leaf flow profile
    leaf_flow_profile!(config, spac, Δt);

    # 2. set up stem flow rate and profile
    stem_flow_profile!(spac, Δt);

    # 3. set up root flow rate and profile
    root_flow_profile!(config, spac, Δt);

    return nothing
);
