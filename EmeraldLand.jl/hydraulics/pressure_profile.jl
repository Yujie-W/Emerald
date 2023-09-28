
#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-May-27: migrate function to version v0.3
#     2022-Jul-12: compute e_crit for leaves
#     2022-Jul-12: compute β for leaves (only for empirical models)
#     2022-Jul-14: update root p_ups from SOIL
#     2022-Oct-20: use add SoilLayer to function variables, because of the removal of SH from RootHydraulics
#     2023-Sep-11: put option update to the SPAC configuration
#     2023-Sep-11: add config to the variable list
#     2023-Sep-11: rename function to critical_flow to xylem_pressure
# To do
#     TODO: add leaf extra-xylary vulnerability curve
#     TODO: compute leaf e_crit
#     TODO: update legacy
#     TODO: make sure to not mixing with top soil that is meant for evaporation
#     TODO: add soil ion concentration
#
#######################################################################################################################################################################################################
