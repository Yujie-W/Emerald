# CLM5 settings
CLM5_PFTG = [0, 2.35, 2.35, 2.35, 4.12, 4.12, 4.45, 4.45, 4.45, 4.7, 4.7, 4.7, 2.22, 5.25, 1.62, 5.79, 5.79] .* sqrt(1000);
CLM5_PFTS = ["not_vegetated",
             "needleleaf_evergreen_temperate",
             "needleleaf_evergreen_boreal",
             "needleleaf_deciduous_boreal",
             "broadleaf_evergreen_tropical",
             "broadleaf_evergreen_temperate",
             "broadleaf_deciduous_tropical",
             "broadleaf_deciduous_temperate",
             "broadleaf_deciduous_boreal",
             "evergreen_shrub",
             "deciduous_temperate_shrub",
             "deciduous_boreal_shrub",
             "c3_arctic_grass",
             "c3_non-arctic_grass",
             "c4_grass",
             "c3_crop",
             "c3_irrigated"];


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: migrate from research repo to Emerald
#
#######################################################################################################################################################################################################
"""

    spac_grids(dts::LandDatasets{FT}) where {FT<:AbstractFloat}

Prepare a matrix of GriddingMachine data to feed SPAC, given
- `dts` `LandDatasets` type data struct

"""
function spac_grids(dts::LandDatasets{FT}) where {FT<:AbstractFloat}
    # read some general data
    _ind_c3 = [2:14;16;17];
    _ccs = read_csv("$(@__DIR__)/../data/CO2-1Y.csv");
    _co2 = _ccs.MEAN[findfirst(_ccs.YEAR .== dts.year)];

    # create a matrix of SPAC
    @tinfo "Preparing a matrix of GriddingMachine data to work on...";
    _mat_gm = Matrix{Union{Nothing,Dict{String,Any}}}(nothing, size(dts.t_lm));
    for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
        if dts.mask_spac[_ilon,_ilat]
            _pfts = dts.t_pft[_ilon,_ilat,:];
            _g = CLM5_PFTG[_ind_c3]' * _pfts[_ind_c3] / sum(_pfts[_ind_c3]);
            _g1 = isnan(_g) ? nanmean(CLM5_PFTG[_ind_c3]) : _g;
            _mat_gm[_ilon,_ilat] = Dict{String,Any}(
                        "CANOPY_HEIGHT" => dts.p_ch[_ilon,_ilat],
                        "CO2"           => _co2,
                        "FT"            => FT,
                        "LATITUDE"      => (_ilat - 0.5) * 180 / size(dts.t_lm,2) - 90,
                        "LMA"           => 1 / dts.p_sla[_ilon,_ilat] / 10,
                        "LONGITUDE"     => (_ilon - 0.5) * 360 / size(dts.t_lm,1) - 180,
                        "MEDLYN_G1"     => _g1,
                        "SOIL_COLOR"    => min(20, max(1, Int(floor(dts.s_cc[_ilon,_ilat])))),
                        "SOIL_N"        => dts.s_n[_ilon,_ilat,:],
                        "SOIL_α"        => dts.s_α[_ilon,_ilat,:],
                        "SOIL_Θr"       => dts.s_Θr[_ilon,_ilat,:],
                        "SOIL_Θr"       => dts.s_Θr[_ilon,_ilat,:],
            );
        end;
    end;

    return _mat_gm
end
