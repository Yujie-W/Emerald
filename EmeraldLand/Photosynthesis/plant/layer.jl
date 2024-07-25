# This file constains function to write the photosynthetic rates into leaf flux struct

#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-25: add leaf_photosynthesis! function for canopy layer
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::Union{GCO₂Mode, PCO₂Mode}; rd_only::Bool = false) where {FT}

Updates leaf photosynthetic rates for the leaf based on leaf stomtal model, given
- `leaf` `Leaf` type structure
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `GCO₂Mode` or `PCO₂Mode` to determine whether to use CO₂ partial pressure or concentration to compute photosynthetic rates
- `rd_only` Whether to compute respiration rate only

"""
function leaf_photosynthesis! end;

# This method takes out stomtal model out and use it to determine whether to apply beta to Vcmax, Jmax, and Rd
leaf_photosynthesis!(leaf::CanopyLayer{FT}, air::AirLayer{FT}, mode::Union{GCO₂Mode, PCO₂Mode}; rd_only::Bool = false) where {FT} = nothing;
