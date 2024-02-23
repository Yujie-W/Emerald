





#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: migrate from research repo to Emerald
#     2023-Jun-15: add non-vegetated land in global simulations
#
#######################################################################################################################################################################################################
"""

    gm_grids(dts::LandDatasets{FT}) where {FT}

Prepare a matrix of GriddingMachine data to feed SPAC, given
- `dts` `LandDatasets` type data struct

"""
function gm_grids(dts::LandDatasets{FT}) where {FT}
    # create a matrix of GriddingMachine data
    # TODO: add a step to verify the input datasets
    @tinfo "Preparing a matrix of GriddingMachine data to work on...";
    mat_gm = Matrix{Union{Nothing,Dict{String,Any}}}(nothing, size(dts.t_lm));
    for ilon in axes(dts.t_lm,1), ilat in axes(dts.t_lm,2)
        if dts.mask_spac[ilon,ilat] || dts.mask_soil[ilon,ilat]
            mat_gm[ilon,ilat] = grid_dict(dts, ilat, ilon);
        end;
    end;

    return mat_gm
end;
