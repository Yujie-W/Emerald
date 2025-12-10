"""

    LATENT_HEAT(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    LATENT_HEAT(spac::BulkSPAC{FT}) where {FT}

Return the latent heat flux per ground area, given
- `config` `SPACConfig` configuration
- `spac` `BulkSPAC` SPAC

"""
function LATENT_HEAT end;

LATENT_HEAT(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = LATENT_HEAT(spac);

LATENT_HEAT(spac::BulkSPAC{FT}) where {FT} = (
    # per ground area
    soils = spac.soils;
    sbulk = spac.soil_bulk;
    soil_le = soils[1].auxil.∂e∂t_le;

    # not area specific
    leaves = spac.plant.leaves;
    leaf_le = FT(0);
    for leaf in leaves
        leaf_le += leaf.energy.auxil.∂e∂t_le;
    end;

    return soil_le + leaf_le / sbulk.trait.area
);


"""

    LONGWAVE_OUT(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    LONGWAVE_OUT(spac::BulkSPAC{FT}) where {FT}

Return the outgoing longwave radiation per ground area, given
- `config` `SPACConfig` configuration
- `spac` `BulkSPAC` SPAC

"""
function LONGWAVE_OUT end;

LONGWAVE_OUT(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = LONGWAVE_OUT(spac);

LONGWAVE_OUT(spac::BulkSPAC{FT}) where {FT} = spac.canopy.structure.auxil.lwꜛ[1];


"""

    NET_LONGWAVE(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    NET_LONGWAVE(spac::BulkSPAC{FT}) where {FT}

Return the net longwave radiation per ground area, given
- `config` `SPACConfig` configuration
- `spac` `BulkSPAC` SPAC

"""
function NET_LONGWAVE end;

NET_LONGWAVE(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = NET_LONGWAVE(spac);

NET_LONGWAVE(spac::BulkSPAC{FT}) where {FT} = (
    # per ground area
    soils = spac.soils;
    sbulk = spac.soil_bulk;
    soil_lw = soils[1].auxil.∂e∂t_lw;

    # not area specific
    leaves = spac.plant.leaves;
    branches = spac.plant.branches;
    canopy_lw = FT(0);
    for leaf in leaves
        canopy_lw += leaf.energy.auxil.∂e∂t_lw;
    end;
    for stem in branches
        canopy_lw += stem.energy.auxil.∂e∂t_lw;
    end;

    return soil_lw + canopy_lw / sbulk.trait.area
);


"""

    NET_SHORTWAVE(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    NET_SHORTWAVE(spac::BulkSPAC{FT}) where {FT}

Return the net shortwave radiation per ground area, given
- `config` `SPACConfig` configuration
- `spac` `BulkSPAC` SPAC

"""
function NET_SHORTWAVE end;

NET_SHORTWAVE(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = NET_SHORTWAVE(spac);

NET_SHORTWAVE(spac::BulkSPAC{FT}) where {FT} = (
    # per ground area
    soils = spac.soils;
    sbulk = spac.soil_bulk;
    soil_sw = soils[1].auxil.∂e∂t_sw;

    # not area specific
    leaves = spac.plant.leaves;
    branches = spac.plant.branches;
    canopy_sw = FT(0);
    for leaf in leaves
        canopy_sw += leaf.energy.auxil.∂e∂t_sw;
    end;
    for stem in branches
        canopy_sw += stem.energy.auxil.∂e∂t_sw;
    end;

    return soil_sw + canopy_sw / sbulk.trait.area
);


"""

    SENSIBLE_HEAT(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    SENSIBLE_HEAT(spac::BulkSPAC{FT}) where {FT}

Return the sensible heat flux per ground area, given
- `config` `SPACConfig` configuration
- `spac` `BulkSPAC` SPAC

"""
function SENSIBLE_HEAT end;

SENSIBLE_HEAT(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = SENSIBLE_HEAT(spac);

SENSIBLE_HEAT(spac::BulkSPAC{FT}) where {FT} = (
    # per ground area
    soils = spac.soils;
    sbulk = spac.soil_bulk;
    soil_sh = soils[1].auxil.∂e∂t_sh;

    # not area specific
    leaves = spac.plant.leaves;
    leaf_sh = FT(0);
    for leaf in leaves
        leaf_sh += leaf.energy.auxil.∂e∂t_sh;
    end;

    return soil_sh + leaf_sh / sbulk.trait.area
);


"""

    SHORTWAVE_OUT(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return the outgoing shortwave radiation per ground area, given
- `config` `SPACConfig` configuration
- `spac` `BulkSPAC` SPAC

"""
function SHORTWAVE_OUT end;

SHORTWAVE_OUT(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = (
    dwl = config.CONSTANTS.SPECTRA.ΔΛ;
    sw_out = spac.canopy.sun_geometry.auxil.e_difꜛ[:,1];

    return (sw_out' * dwl) / 1000
);


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
