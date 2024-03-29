#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-11: add function to compute PAR above canopy
#
#######################################################################################################################################################################################################
"""

    PAR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Return the PAR above canopy per ground area, given
- `config` `SPACConfiguration` SPAC configuration
- `spac` `BulkSPAC` SPAC

"""
function PAR end;

PAR(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = (
    (; SPECTRA) = config;
    rad_sw = spac.meteo.rad_sw;

    ppfd = photon.(SPECTRA.Λ_PAR, (rad_sw.e_dif + rad_sw.e_dir)[SPECTRA.IΛ_PAR]) .* 1000;

    return ppfd' * SPECTRA.ΔΛ_PAR
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Oct-19: add function to compute canopy integrated PPAR
#     2023-May-19: use δlai per canopy layer
#
#######################################################################################################################################################################################################
"""

    PPAR(spac::BulkSPAC{FT}) where {FT}

Return the canopy integrated PPAR per ground area, given
- `spac` `BulkSPAC` SPAC

"""
function PPAR end;

PPAR(spac::BulkSPAC{FT}) where {FT} = (
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    # compute GPP
    ppar::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        ppar += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.ppar_sunlit) +
                (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.ppar_shaded) * canopy.structure.state.δlai[i];
    end;

    return ppar
);
