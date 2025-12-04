"""
Struct that contains the trait and state variables.
"""
Base.@kwdef mutable struct LeafTraitsStates{FT<:AbstractFloat}
    # leaf pigments
    "Anthocyanin content `[μg cm⁻²]`"
    ant::FT = 0
    "Senescent material (brown pigments) fraction `[-]`"
    brown::FT = 0
    "Chlorophyll a and b content `[μg cm⁻²]`"
    cab::FT = 40
    "Carotenoid content `[μg cm⁻²]`"
    car::FT = 40 / 7
    "Carbon-based constituents in lma `[g cm⁻²]`"
    cbc::FT = 0
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    lma::FT = 0.012
    "Leaf mesophyll structural parameter that describes the number of thin layers with a leaf"
    meso_n::FT = 1.4
    "Protein content in lma (pro = lma - cbc) `[g cm⁻²]`"
    pro::FT = 0

    # leaf pigments related states
    "Zeaxanthin fraction in Carotenoid (1=all Zeaxanthin, 0=all Violaxanthin) `[-]`"
    f_zeax::FT = 0
    "Fraction of carotenoid absorption into SIF `[-]`"
    ϕ_car::FT = 1
    "Fraction of carotenoid absorption into PPAR `[-]`"
    ϕ_car_ppar::FT = 1

    # leaf width
    "leaf width `[m]`"
    width::FT = 0.05

    # longwave radiation
    "Broadband thermal reflectance, related to blackbody emittance `[-]`"
    ρ_lw::FT = 0.01
    "Broadband thermal transmission, related to blackbody emittance `[-]`"
    τ_lw::FT = 0.01
end;
