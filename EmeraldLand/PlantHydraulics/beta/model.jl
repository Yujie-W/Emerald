#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-30: migrate function from older version StomataModels.jl
#     2022-Jul-01: add method to tune stomatal opening based on relative hydraulic conductance at leaf xylem end;
#     2022-Jul-01: add method to tune stomatal opening based on relative hydraulic conductance of the soil
#     2022-Jul-01: add method to tune stomatal opening based on soil potential or leaf pressure
#     2022-Nov-18: force the beta to be within (0,1]
#
#######################################################################################################################################################################################################
"""

    β_factor(f::Function, vc::AbstractXylemVC{FT}, x_25::FT) where {FT}
    β_factor(f::Function, vc::AbstractSoilVC{FT}, x_25::FT) where {FT}
    β_factor(f::Function, x_25::FT) where {FT}

Return the β factor based on relative conductance or soil potential/pressure, given
- `f` Function to translate relative k to β, for example f(x) = x, f(x) = x², and f(x) = sqrt(x) for x in [0,1]
- `vc` Leaf vulnerability curve or soil vulnerability curve (moisture retention curve)
- `x_25` Leaf xylem pressure corrected to 25 °C, soil water potential corrected to 25 °C (forcing on roots, note that this function may not be useful for plants with salt stress), or soil water
    content.

Here, we pass a `f` (Function) into the function call so that we can easily customized the beta function without hacking the source code.

"""
function β_factor end;

β_factor(f::Function, vc::AbstractXylemVC{FT}, x_25::FT) where {FT} = max(eps(FT), min(1, f(relative_soil_k(vc, x_25))));

β_factor(f::Function, vc::AbstractSoilVC{FT}, x_25::FT) where {FT} = max(eps(FT), min(1, f(relative_soil_k(vc, true, x_25))));

β_factor(f::Function, x_25::FT) where {FT} = max(eps(FT), min(1, f(x_25)));
