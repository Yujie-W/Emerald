
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

=#
