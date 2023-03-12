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

    spac_grids(dts::LandDatasets{FT}; threads::Int = 12) where {FT<:AbstractFloat}

Prepare a matrix of SPAC, given
- `dts` `LandDatasets` type data struct
- `threads` Number of threadings to use

"""
function spac_grids(dts::LandDatasets{FT}; threads::Int = 12) where {FT<:AbstractFloat}
    add_threads!(threads);

    # read some general data
    _ind_c3 = [2:14;16;17];
    _ccs = read_csv("$(@__DIR__)/../data/CO2-1Y.csv");
    _co2 = _ccs.MEAN[findfirst(_ccs.YEAR .== dts.year)];

    # create a matrix of SPAC
    @tinfo "Preparing a matrix of SPAC to work on...";
    _mat_spac = Matrix{Union{Nothing,MonoMLTreeSPAC}}(nothing, size(dts.t_lm));
    _params = [];
    for _ilon in axes(dts.t_lm,1), _ilat in axes(dts.t_lm,2)
        if dts.mask_spac[_ilon,_ilat]
            push!(_params, [_ilon,_ilat]);
        end;
    end;

    @inline linear_p_soil(x) = min(1, max(eps(FT), 1 + x / 5));
    @inline create_spac(param) = (
        (_ilon,_ilat) = param;
        # create a SPAC to work on
        _z_canopy = max(FT(0.1), dts.p_ch[_ilon,_ilat]);
        _spac = MonoMLTreeSPAC{FT}(
                    DIM_AIR      = 25,
                    DIM_LAYER    = 10,
                    DIM_ROOT     = 4,
                    LATITUDE     = (_ilat - 0.5) * 180 / size(dts.t_lm,2),
                    LONGITUDE    = (_ilon - 0.5) * 360 / size(dts.t_lm,1),
                    LEAVES_INDEX = collect(11:20),
                    ROOTS_INDEX  = collect(1:4),
                    Z            = [-2, _z_canopy/2, _z_canopy],
                    Z_AIR        = collect(0:21) * _z_canopy / 20,
                    SOIL         = Soil{FT}(DIM_SOIL = 4, COLOR = min(20, max(1, Int(floor(dts.s_cc[_ilon,_ilat])))), ZS = [0, -0.1, -0.35, -1, -3]));

        # update soil type information per layer
        for _i in eachindex(_spac.SOIL.LAYERS)
            # TODO: add a line to parameterize K_MAX
            _spac.SOIL.LAYERS[_i].VC.α = dts.s_α[_ilon,_ilat,_i];
            _spac.SOIL.LAYERS[_i].VC.N = dts.s_n[_ilon,_ilat,_i];
            _spac.SOIL.LAYERS[_i].VC.M = 1 - 1 / _spac.SOIL.LAYERS[_i].VC.N;
            _spac.SOIL.LAYERS[_i].VC.Θ_RES = dts.s_Θr[_ilon,_ilat,_i];
            _spac.SOIL.LAYERS[_i].VC.Θ_SAT = dts.s_Θs[_ilon,_ilat,_i];
        end;

        # set hydraulic traits to very high so as to not triggering NaN (they do not impact result anyway)
        for _organ in [_spac.LEAVES; _spac.BRANCHES; _spac.TRUNK; _spac.ROOTS]
            _organ.HS.VC.B = 3;
            _organ.HS.VC.C = 1;
        end;

        # update leaf mass per area and stomtal model
        _pfts = dts.t_pft[_ilon,_ilat,:];
        _g = CLM5_PFTG[_ind_c3]' * _pfts[_ind_c3] / sum(_pfts[_ind_c3]);
        _g1 = isnan(_g) ? nanmean(CLM5_PFTG[_ind_c3]) : _g;
        _bt = BetaFunction{FT}(FUNC = linear_p_soil, PARAM_X = BetaParameterPsoil(), PARAM_Y = BetaParameterG1());
        for _leaves in _spac.LEAVES
            _leaves.BIO.lma = 1 / dts.p_sla[_ilon,_ilat] / 10;
            _leaves.SM = MedlynSM{FT}(G0 = 0.005, G1 = _g1, β = _bt);
        end;

        # update the vcmax for C3 model
        update!(_spac; vcmax = dts.p_vcm[_ilon,_ilat,1], vcmax_expo = 0.3);

        # sync the environmental conditions per layer for CO₂ concentration
        for _alayer in _spac.AIR
            update!(_alayer; f_CO₂ = _co2);
        end;

        # initialize the spac
        initialize!(_spac);

        return (_ilon,_ilat,_spac)
    );

    _thread_spacs = @showprogress pmap(create_spac, _params);
    for _thread_spac in _thread_spacs
        _ilon,_ilat,_spac = _thread_spac;
        _mat_spac[_ilon,_ilat] = _spac;
    end;

    return _mat_spac
end
