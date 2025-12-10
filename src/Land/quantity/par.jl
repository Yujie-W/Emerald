"""

    APAR(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}
    APAR(spac::BulkSPAC{FT}) where {FT}

Return the canopy integrated APAR per ground area, given
- `config` `SPACConfig` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function APAR end;

APAR(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = APAR(spac);

APAR(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # compute GPP
    apar::FT = 0;
    for irt in eachindex(leaves)
        ilf = n_layer + 1 - irt;
        apar += leaves[ilf].flux.auxil.apar' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt) * canopy.structure.trait.δlai[irt];
    end;

    return apar
);


"""

    PAR(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT}

Return the PAR above canopy per ground area, given
- `config` `SPACConfig` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function PAR end;

PAR(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = (
    (; SPECTRA) = config.CONSTANTS;
    rad_sw = spac.meteo.rad_sw;

    ppfd::FT = 0;
    for i_par in eachindex(SPECTRA.Λ_PAR)
        ppfd += energy_to_photon(SPECTRA.Λ_PAR[i_par], rad_sw.e_dif[SPECTRA.IΛ_PAR[i_par]] + rad_sw.e_dir[SPECTRA.IΛ_PAR[i_par]]) * 1000 * SPECTRA.ΔΛ_PAR[i_par];
    end;

    return ppfd
);


"""

    PPAR(spac::BulkSPAC{FT}) where {FT}

Return the canopy integrated PPAR per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function PPAR end;

PPAR(config::SPACConfig{FT}, spac::BulkSPAC{FT}) where {FT} = PPAR(spac);

PPAR(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    n_layer = length(leaves);

    # compute GPP
    ppar::FT = 0;
    for irt in eachindex(leaves)
        ilf = n_layer + 1 - irt;
        ppar += leaves[ilf].flux.auxil.ppar' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt) * canopy.structure.trait.δlai[irt];
    end;

    return ppar
);
