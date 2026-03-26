"""

    T_SKIN(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    T_SKIN(spac::BulkSPAC{FT}) where {FT}

Return the skin temperature, given
- `config` `SPACConfig` configuration
- `spac` `BulkSPAC` SPAC

"""
function T_SKIN end;

T_SKIN(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = T_SKIN(spac);

T_SKIN(spac::BulkSPAC{FT}) where {FT} = (
    lw_out = LONGWAVE_OUT(spac);

    return FT( (lw_out / FT(0.98) / K_STEFAN()) ^ (1/4) )
);


"""

    T_SOIL(config::SPACConfig{FT}, spac::BulkSPAC{FT}, ilayer::Int) where {FT}
    T_SOIL(spac::BulkSPAC{FT}, ilayer::Int) where {FT}

Return the soil temperature of a given soil layer (counting from the top), given
- `config` `SPACConfig` configuration
- `spac` `BulkSPAC` SPAC
- `ilayer` soil layer index

"""
function T_SOIL end;

T_SOIL(config::SPACConfig{FT}, spac::BulkSPAC{FT}, ilayer::Int) where {FT} = T_SOIL(spac, ilayer);

T_SOIL(spac::BulkSPAC{FT}, ilayer::Int) where {FT} = (
    return spac.soils[ilayer].auxil.t
);
