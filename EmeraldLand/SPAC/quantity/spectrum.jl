
#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-13: add method to interpolate the spectrum
#     2022-Jun-13: add method to interpolate the spectrum via multiple steps
#     2024-Nov-13: move the method of read_spectrum for single x as interpolate_data
#
#######################################################################################################################################################################################################
"""

    read_spectrum(x::Vector{FT}, y::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT}
    read_spectrum(x::Vector{FT}, y::Vector{FT}, w::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT}

Return the spectrum value at target wavelength bin, given
- `x` X-axis of the spectrum
- `y` Y-axis of the spectrum
- `w` Weight of the spectrum
- `x₁` Lower x boundary
- `x₂` Upper x boundary
- `steps` The incremental Δx is `(x₂ - x₁) / steps`

"""
function read_spectrum end;

read_spectrum(x::Vector{FT}, y::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT} = (
    ys = 0;
    δx = (x₂ - x₁) / steps;
    for i in 1:(steps+1)
        xi = x₁ + (i - 1) * δx;
        ys += interpolate_data(x, y, xi);
    end;

    return ys / (steps + 1)
);

read_spectrum(x::Vector{FT}, y::Vector{FT}, w::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT} = (
    ys = 0;
    ws = 0;
    δx = (x₂ - x₁) / steps;
    for i in 1:(steps+1)
        xi = x₁ + (i - 1) * δx;
        ys += interpolate_data(x, y, xi);
        ws += interpolate_data(x, w, xi);
    end;

    return ys / ws
);
