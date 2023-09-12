#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-30: migrate function from older version StomataModels.jl
#     2022-Jul-01: add method to tune stomatal opening based on relative hydraulic conductance at leaf xylem end
#     2022-Jul-01: add method to tune stomatal opening based on relative hydraulic conductance of the soil
#     2022-Jul-01: add method to tune stomatal opening based on soil potential or leaf pressure
#     2022-Jul-12: move function from StomataModels.jl to PlantHydraulics.jl
#     2022-Nov-18: force the beta to be within (0,1]
#     2023-Aug-27: add nan check for beta calculation of empirical models (optimality models have a beta = nan)
#
#######################################################################################################################################################################################################
"""

    β_factor(f::Function, vc::AbstractXylemVC{FT}, x_25::FT) where {FT<:AbstractFloat}
    β_factor(f::Function, vc::AbstractSoilVC{FT}, x_25::FT) where {FT<:AbstractFloat}
    β_factor(f::Function, x_25::FT) where {FT<:AbstractFloat}

Return the β factor based on relative conductance or soil potential/pressure, given
- `f` Function to translate relative k to β, for example f(x) = x, f(x) = x², and f(x) = sqrt(x) for x in [0,1]
- `vc` Leaf vulnerability curve or soil vulnerability curve (moisture retention curve)
- `x_25` Leaf xylem pressure corrected to 25 °C, soil water potential corrected to 25 °C (forcing on roots, note that this function may not be useful for plants with salt stress), or soil water
    content.

"""
function β_factor end

β_factor(f::Function, vc::AbstractXylemVC{FT}, x_25::FT) where {FT<:AbstractFloat} = (
    _β = FT(max(eps(FT), min(1, f(relative_hydraulic_conductance(vc, x_25)))));

    if isnan(_β)
        @info "Debugging" vc x_25;

        return error("Computed β is NaN")
    else
        return _β
    end
);

β_factor(f::Function, vc::AbstractSoilVC{FT}, x_25::FT) where {FT<:AbstractFloat} = (
    _β = FT(max(eps(FT), min(1, f(relative_hydraulic_conductance(vc, true, x_25)))));

    if isnan(_β)
        @info "Debugging" vc x_25;

        return error("Computed β is NaN")
    else
        return _β
    end
);

β_factor(f::Function, x_25::FT) where {FT<:AbstractFloat} = (
    _β = FT(max(eps(FT), min(1, f(x_25))));

    if isnan(_β)
        @info "Debugging" x_25;

        return error("Computed β is NaN")
    else
        return _β
    end
);
