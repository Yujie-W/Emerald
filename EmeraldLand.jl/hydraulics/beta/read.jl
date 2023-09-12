#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-11: rename the methods from β_factor to read_β
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
function read_β end

read_β(sm::AbstractStomataModel{FT}) where {FT<:AbstractFloat} = FT(NaN);

read_β(sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}) where {FT<:AbstractFloat} = sm.β.β₁;
