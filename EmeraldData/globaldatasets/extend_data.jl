#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: migrate from research repo to Emerald
#     2023-Mar-11: add method to extend the LandDatasets
#     2023-Aug-25: move method on vectors to EmeraldMath.jl
#
#######################################################################################################################################################################################################
"""

    extend_data!(dts::LandDatasets{FT}) where {FT}

Gap fill the data linearly, given
- `dts` LandDatasets struct

"""
function extend_data!(dts::LandDatasets{FT}) where {FT}
    # determine where to fill based on land mask and lai
    for ilon in axes(dts.t_lm,1), ilat in axes(dts.t_lm,2)
        if (dts.t_lm[ilon,ilat] > 0) && (nanmax(dts.p_lai[ilon,ilat,:]) > 0)
            dts.mask_spac[ilon,ilat] = true;
            mask_lai = isnan.(dts.p_lai[ilon,ilat,:]);
            dts.p_lai[ilon,ilat,mask_lai] .= 0;
        elseif (dts.t_lm[ilon,ilat] > 0)
            dts.mask_soil[ilon,ilat] = true;
            mask_lai = isnan.(dts.p_lai[ilon,ilat,:]);
            dts.p_lai[ilon,ilat,mask_lai] .= 0;
        end;
    end;

    # iterate the fieldnames
    for fn in fieldnames(typeof(dts))
        if !(fn in [:p_lai, :t_ele, :t_lm, :t_pft, :mask_soil, :mask_spac])
            data = getfield(dts, fn);
            if data isa Array
                # extend the data first based on interpolations
                for ilon in axes(dts.t_lm,1), ilat in axes(dts.t_lm,2)
                    tmp = data[ilon,ilat,:];
                    interpolate_data!(tmp);
                    data[ilon,ilat,:] .= tmp;
                end;

                # fill the NaNs with nanmean of the rest
                mask_mean = dts.mask_spac .&& isnan.(data);
                data[mask_mean] .= nanmean(data);
            end;
        end;
    end;

    return nothing
end;
