# This file contains function to compute the sun geometry related parameters (apart from sensor geometry)

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Oct-10: add function sun_geometry! (run per solar zenith angle)
#
#######################################################################################################################################################################################################
"""

    sun_geometry!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT}

Update sun geometry related auxiliary variables, given
- `config` SPAC configuration
- `can` SPAC canopy

"""
function sun_geometry!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}) where {FT}
    (; DIM_LAYER, Θ_AZI, Θ_INCL) = config;

    # extinction coefficients for the solar radiation
    for i in eachindex(Θ_INCL)
        Cs = cosd(Θ_INCL[i]) * cosd(can.sun_geometry.state.sza);
        Ss = sind(Θ_INCL[i]) * sind(can.sun_geometry.state.sza);
        βs = (Cs >= Ss ? FT(π) : acos(-Cs/Ss));
        can.sun_geometry.auxil.Cs_incl[i] = Cs;
        can.sun_geometry.auxil.Ss_incl[i] = Ss;
        can.sun_geometry.auxil.βs_incl[i] = βs;
        can.sun_geometry.auxil.ks_incl[i] = 2 / FT(π) / cosd(can.sun_geometry.state.sza) * (Cs * (βs - FT(π)/2) + Ss * sin(βs));
    end;
    can.sun_geometry.auxil.ks = can.structure.state.p_incl' * can.sun_geometry.auxil.ks_incl;

    # compute the scattering weights for diffuse/direct -> diffuse for backward and forward scattering
    can.sun_geometry.auxil.ddb = (1 + can.structure.auxil.bf) / 2;
    can.sun_geometry.auxil.ddf = (1 - can.structure.auxil.bf) / 2;
    can.sun_geometry.auxil.sdb = (can.sun_geometry.auxil.ks + can.structure.auxil.bf) / 2;
    can.sun_geometry.auxil.sdf = (can.sun_geometry.auxil.ks - can.structure.auxil.bf) / 2;

    # compute the fs and fs_abs matrices
    for i in eachindex(Θ_AZI)
        view(can.sun_geometry.auxil.fs,:,i) .= can.sun_geometry.auxil.Cs_incl .+ can.sun_geometry.auxil.Ss_incl .* cosd(Θ_AZI[i]);
    end;
    can.sun_geometry.auxil.fs ./= cosd(can.sun_geometry.state.sza);
    can.sun_geometry.auxil.fs_abs .= abs.(can.sun_geometry.auxil.fs);

    # compute the sunlit leaf fraction
    can.sun_geometry.auxil.ps = exp.(can.sun_geometry.auxil.ks .* can.structure.auxil.ci * can.structure.state.lai .* can.structure.auxil.x_bnds);
    for i in 1:DIM_LAYER
        ilai = can.structure.auxil.ci * can.structure.state.δlai[i];
        kscilai = can.sun_geometry.auxil.ks * can.structure.auxil.ci * can.structure.state.lai;
        can.sun_geometry.auxil.p_sunlit[i] = 1 / (can.sun_geometry.auxil.ks * ilai) * (exp(kscilai * can.structure.auxil.x_bnds[i]) - exp(kscilai * can.structure.auxil.x_bnds[i+1]));
    end;

    return nothing
end;
