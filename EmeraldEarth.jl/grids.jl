
function spac_grids(dts::LandDatasets{FT}) where {FT<:AbstractFloat}
    _mat_spac = Matrix{Union{Nothing,MonoMLTreeSPAC}}(undef, size(dts.t_lm));

    return _mat_spac
end
