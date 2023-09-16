#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-15: add function to compute effective interface reflectance from known layer reflectance, transmittance, and absorption (the two functions need to be used at the same time)
#
# Algorithm
#     here we want to solve the following equation for ρ_12 and ρ_21 (decouple the two)
#         τ_eff =                 (1 - ρ_12) * τ_all * (1 - ρ_21) / (1 - τ_all ^ 2 * ρ_21 ^ 2)
#         ρ_eff = x + τ_all * y * (1 - ρ_12) * τ_all * (1 - ρ_21) / (1 - τ_all ^ 2 * ρ_21 ^ 2)
#     let x = ρ_12
#         y = ρ_21
#         t = τ_eff
#         r = ρ_eff
#         k = τ_all
#     then we have
#         t =             (1 - x) * k * (1 - y) / (1 - k ^ 2 * y ^ 2)
#         r = x + k * y * (1 - x) * k * (1 - y) / (1 - k ^ 2 * y ^ 2)
#     solve x and y using wolframalpha.com and we have
#         x = (t * (k - t) + r^2 - r) / (k * t + r  - 1)
#         y = (k * (r - 1) + t) / k / (k * t + r - 1)
#
#######################################################################################################################################################################################################
"""

    effective_ρ_12(ρ_eff::FT, τ_eff::FT, τ_all::FT) where {FT}

Return the effective air-water interface reflectance, given
- `ρ_eff` effective reflectance of the layer
- `τ_eff` effective transmittance of the layer
- `τ_all` transmittance within the layer

"""
function effective_ρ_12(ρ_eff::FT, τ_eff::FT, τ_all::FT) where {FT}
    return (τ_eff * (τ_all - τ_eff) + ρ_eff^2 - ρ_eff) / (τ_all * τ_eff + ρ_eff - 1)
end;


"""

    effective_ρ_21(ρ_eff::FT, τ_eff::FT, τ_all::FT) where {FT}

Return the effective water-air interface reflectance, given
- `ρ_eff` effective reflectance of the layer
- `τ_eff` effective transmittance of the layer
- `τ_all` transmittance within the layer

"""
function effective_ρ_21(ρ_eff::FT, τ_eff::FT, τ_all::FT) where {FT}
    return (τ_all * (ρ_eff - 1) + τ_eff) / τ_all / (τ_all * τ_eff + ρ_eff - 1)
end;
