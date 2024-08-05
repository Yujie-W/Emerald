
#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-13: add method to interpolate the spectrum
#     2022-Jun-13: add method to interpolate the spectrum via multiple steps
#
#######################################################################################################################################################################################################
"""

    read_spectrum(x::Vector{FT}, y::Vector{FT}, target::FT) where {FT}
    read_spectrum(x::Vector{FT}, y::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT}
    read_spectrum(x::Vector{FT}, y::Vector{FT}, w::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT}

Return the spectrum value at target wavelength bin, given
- `x` X-axis of the spectrum
- `y` Y-axis of the spectrum
- `w` Weight of the spectrum
- `target` Target x value
- `x₁` Lower x boundary
- `x₂` Upper x boundary
- `steps` The incremental Δx is `(x₂ - x₁) / steps`

"""
function read_spectrum end;

read_spectrum(x::Vector{FT}, y::Vector{FT}, target::FT) where {FT} = (
    @assert length(x) == length(y) "Dimensions of provided spectrum x and y must match!";
    @assert x[1] <= target <= x[end] "Target wavelength must be within the range provided spectum!";

    # iterate through the spectrum and find the index
    ind = 0;
    for i in 1:length(x)-1
        if x[i] <= target <= x[i+1]
            ind = i;
            break;
        end;
    end;

    return ((x[ind+1] - target) * y[ind] + (target - x[ind]) * y[ind+1]) / (x[ind+1] - x[ind])
);

read_spectrum(x::Vector{FT}, y::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT} = (
    ys = 0;
    δx = (x₂ - x₁) / steps;
    for i in 1:(steps+1)
        xi = x₁ + (i - 1) * δx;
        ys += read_spectrum(x, y, xi);
    end;

    return ys / (steps + 1)
);

read_spectrum(x::Vector{FT}, y::Vector{FT}, w::Vector{FT}, x₁::FT, x₂::FT; steps::Int = 2) where {FT} = (
    ys = 0;
    ws = 0;
    δx = (x₂ - x₁) / steps;
    for i in 1:(steps+1)
        xi = x₁ + (i - 1) * δx;
        ys += read_spectrum(x, y, xi);
        ws += read_spectrum(x, w, xi);
    end;

    return ys / ws
);
