#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: Move the structure from Photosynthesis.jl, only P_A and P_O2 for now
#     2022-Jan-14: rename P_A to P_AIR, P_O2 to P_O₂
#     2022-Jan-14: add p_CO₂ to the structure
#     2022-Mar-09: add t, p_H₂O, p_H₂O_sat, rh, and wind fields
#     2022-Jul-13: remove fields p_H₂O_sat and rh to avoid update issues
#     2022-Jul-20: remove fields P_O₂ to avoid update issues
#     2022-Jul-20: add fields: Z, ΔZ, e, n_CO₂, n_H₂O, ∂e∂t, ∂CO₂∂t, and ∂H₂O∂t
#     2023-Mar-11: add field f_CO₂ for CO₂ concentration
#     2023-Jun-13: add fields for CH₄, N₂, O₂ partial pressures and moles
#     2023-Jul-06: correct the values for CH₄, N₂, O₂ moles
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores air layer information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct AirLayer{FT<:AbstractFloat}
    # Location and geometry of the air layer
    "Mean height of the layer `[m]`"
    Z::FT = 0.5
    "Layer thickness `[m]`"
    ΔZ::FT = 1

    # Parameters that are not supposed to change with time
    "Atmospheric pressure `[Pa]`"
    P_AIR::FT = P_ATM(FT)

    # Prognostic variables (not used for ∂y∂t)
    "CO₂ concentration `[ppm]`"
    f_CO₂::FT = 400
    "CH₄ partial pressure `[Pa]`"
    p_CH₄::FT = 0
    "CO₂ partial pressure `[Pa]`"
    p_CO₂::FT = P_AIR * f_CO₂ * 1e-6
    "H₂O partial pressure `[Pa]`"
    p_H₂O::FT = 1500
    "N₂ partial pressure `[Pa]`"
    p_N₂::FT = P_AIR * 0.79
    "O₂ partial pressure `[Pa]`"
    p_O₂::FT = P_AIR * 0.209
    "Temperature `[K]`"
    t::FT = T₂₅(FT)
    "Wind speed `[m s⁻¹]`"
    wind::FT = 1

    # Prognostic variables (used for ∂y∂t)
    "Total energy within the air layer `[J m⁻²]`"
    Σe::FT = CP_D_MOL(FT) * (P_AIR - p_H₂O) * ΔZ / GAS_R(FT) + CP_V_MOL(FT) * p_H₂O * ΔZ / GAS_R(FT)
    "Mole of CH₄ per surface area `[mol m⁻²]`"
    n_CH₄::FT = p_CH₄ * ΔZ / (GAS_R(FT) * t)
    "Mole of CO₂ per surface area `[mol m⁻²]`"
    n_CO₂::FT = p_CO₂ * ΔZ / (GAS_R(FT) * t)
    "Mole of H₂O per surface area `[mol m⁻²]`"
    n_H₂O::FT = p_H₂O * ΔZ / (GAS_R(FT) * t)
    "Mole of N₂ per surface area `[mol m⁻²]`"
    n_N₂::FT = p_N₂ * ΔZ / (GAS_R(FT) * t)
    "Mole of O₂ per surface area `[mol m⁻²]`"
    n_O₂::FT = p_O₂ * ΔZ / (GAS_R(FT) * t)
    "Marginal increase in total energy `[J m⁻² s⁻¹]`"
    ∂e∂t::FT = 0
    "Marginal increase in total moles of CO₂ `[mol m⁻² s⁻¹]`"
    ∂CO₂∂t::FT = 0
    "Marginal increase in total moles of H₂O `[mol m⁻² s⁻¹]`"
    ∂H₂O∂t::FT = 0
end
