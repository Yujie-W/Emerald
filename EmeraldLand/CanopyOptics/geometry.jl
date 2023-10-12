
#=
canopy_optical_properties!(config::SPACConfiguration{FT}, can::MultiLayerCanopy{FT}, leaves::Vector{Leaf{FT}}, sbulk::SoilBulk{FT}) where {FT} = (
    (; DIM_LAYER) = config;
    (; OPTICS) = can;
    @assert length(leaves) == DIM_LAYER "Number of leaves must be equal to the canopy layers!";

    if can.structure.state.lai == 0
        OPTICS.ρ_dd  .= 0;
        OPTICS.ρ_lw  .= 0;
        OPTICS.ρ_sd  .= 0;
        OPTICS.τ_dd  .= 0;
        OPTICS.τ_lw  .= 0;
        OPTICS.τ_sd  .= 0;
        OPTICS._τ_ss .= 0;
        OPTICS.ρ_dd[:,end] .= sbulk.auxil.ρ_sw;
        OPTICS.ρ_sd[:,end] .= sbulk.auxil.ρ_sw;
        OPTICS.ρ_lw[end] = sbulk.auxil.ρ_lw;

        return nothing
    end;

    return nothing
);

shortwave radiation

    if can.structure.state.lai == 0
        for i in 1:(DIM_LAYER+1)
            can.sun_geometry.auxil.e_dirꜜ[:,i] .= rad.e_dir;
            can.sun_geometry.auxil.e_difꜜ[:,i] .= rad.e_dif;
            can.sun_geometry.auxil.e_difꜛ[:,i] .= 0;
            RADIATION.e_v[:,i] .= 0;
        end;
        can.sun_geometry.auxil.e_difꜛ[:,end] .= view(OPTICS.ρ_sd,:,DIM_LAYER+1) .* view(can.sun_geometry.auxil.e_dirꜜ,:,DIM_LAYER+1) .+ view(OPTICS.ρ_dd,:,DIM_LAYER+1) .* view(can.sun_geometry.auxil.e_difꜜ,:,DIM_LAYER+1);
        RADIATION.e_v[:,end] .= view(can.sun_geometry.auxil.e_difꜛ,:,DIM_LAYER+1);
        RADIATION.e_o .= view(can.sun_geometry.auxil.e_difꜛ,:,DIM_LAYER+1) ./ FT(π);
        RADIATION.albedo .= RADIATION.e_o * FT(π) ./ (rad.e_dir .+ rad.e_dif);

        RADIATION._par_shaded .= photon.(SPECTRA.Λ_PAR, view(rad.e_dif,SPECTRA.IΛ_PAR)) .* 1000;
        RADIATION._par_sunlit .= photon.(SPECTRA.Λ_PAR, view(rad.e_dir,SPECTRA.IΛ_PAR)) .* 1000;
        RADIATION.par_in_diffuse = RADIATION._par_shaded' * SPECTRA.ΔΛ_PAR;
        RADIATION.par_in_direct = RADIATION._par_sunlit' * SPECTRA.ΔΛ_PAR;
        RADIATION.par_in = RADIATION.par_in_diffuse + RADIATION.par_in_direct;

        return nothing
    end;

longwave radiation

    if can.structure.state.lai == 0
        _r_lw_soil = K_STEFAN(FT) * (1 - sbulk.auxil.ρ_lw) * soil.auxil.t ^ 4;
        RADIATION.r_lw .= 0;
        RADIATION.r_net_lw .= 0;
        RADIATION.r_lw_up .= rad * sbulk.auxil.ρ_lw + _r_lw_soil;
        sbulk.auxil.r_net_lw = rad * (1 - sbulk.auxil.ρ_lw) - _r_lw_soil;

        return nothing
    end;

=#
