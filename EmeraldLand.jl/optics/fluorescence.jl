#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-16: add function to update the SIF conversion matrix of the first layer
#     2023-Sep-16: clear vec_b and vec_f before the calculation (for the case of recalculating the SIF conversion matrix)
#
#######################################################################################################################################################################################################
"""

    layer_1_sif_vec!(τ_i_θ::FT, τ_i_12::FT, τ_i_21::FT, τ_sub::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, f_sife::FT, τ_sub_sif::SubArray, vec_b::SubArray, vec_f::SubArray, ϕ::SubArray, N::Int) where {FT}

Update the SIF conversion matrix of the first layer, given
- `τ_i_θ` transmittance of the incoming radiation at the air-water interface
- `τ_i_12` transmittance of the isotropic radiation at the air-water interface
- `τ_i_21` transmittance of the isotropic radiation at the water-air interface
- `τ_sub` transmittance within a sublayer
- `τ_θ` transmittance of the incoming radiation across the leaf layer
- `ρ_1` reflectance of the first layer (1)
- `ρ_2` reflectance of the second layer (n-1)
- `f_sife` SIF excitation scaling factor (f_cab + f_car * ϕ_car)
- `τ_sub_sif` transmittance within a sublayer for SIF
- `vec_b` SIF vector backwards (SubArray of a Matrix)
- `vec_f` SIF vector forwards (SubArray of a Matrix)
- `ϕ` SIF emission PDF (SubArray of a Vector)
- `N` number of sublayers of each layer

"""
function layer_1_sif_vec!(τ_i_θ::FT, τ_i_12::FT, τ_i_21::FT, τ_sub::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, f_sife::FT, τ_sub_sif::SubArray, vec_b::SubArray, vec_f::SubArray, ϕ::SubArray, N::Int) where {FT}
    # parameters required for the calculation that can be derived from the input parameters
    ρ_i_21 = 1 - τ_i_21;
    α_sub = 1 - τ_sub;
    τ_all = τ_sub ^ N;

    # clear the SIF vector
    vec_b .= 0;
    vec_f .= 0;

    # 1. note that we do redo the calculation of the reflectance and transmittance within the leaf layer so that we can better model the sif emission and escape
    #    here the radiation is directly from the environment, as the reflection from the first layer will not hit back, no scaling of the light source is required
    #    however, the rescaling from the reflectance within the layer is required (we will do it at the last step)
    #    for the emitted SIF, we need to account for the reabsorption within the leaf layer
    #        - the SIF up will be downscaled by _t_sub ^ (_i - 0.5) times
    #        - the SIF down will be downscaled by _t_sub ^ (N - _i + 0.5) times
    #    ϕ_sif to account for the contribution to SIF excitation
    #    now the radiation goes from up to down
    rad_i = τ_i_θ * f_sife;
    for i in 1:N
        vec_b .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 2. then the radiation goes from down to up after hitting the water-air interface
    rad_i *= ρ_i_21;
    for i in N:-1:1
        vec_b .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 3. we also need to acount for the radiation reflected by the lower n-1 layer
    #    as this radiation is reflected back and forth, here we use the total radiation that goes in the direction into the first layer
    #        τ_θ * ρ_2 / (1 - ρ_1 * ρ_2)
    #    after accounting for the reflection and transmission within the layer
    #        the total radiation that transmits through is τ_θ * ρ_2 / (1 - ρ_1 * ρ_2) * τ_1
    #        the total radiation that reflects back is τ_θ * ρ_2 / (1 - ρ_1 * ρ_2) * ρ_1
    #    in this case, the total radiation is conserved
    #    now the radiation goes from down to up
    rad_i = τ_θ * ρ_2 / (1 - ρ_1 * ρ_2) * τ_i_12 * f_sife;
    for i in N:-1:1
        vec_b .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 4. then the radiation goes from up to down
    rad_i *= ρ_i_21;
    for i in 1:N
        vec_b .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 5. rescale it by 1 / (1 - ρ_21 * t_all * ρ_21 * t_all) for the SIFE
    denom = 1 - τ_all * ρ_i_21 * τ_all * ρ_i_21;
    vec_b ./= denom;
    vec_f ./= denom;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-16: add function to update the SIF conversion maxtri of the n-1 layer
#
#######################################################################################################################################################################################################
"""

    layer_2_sif_vec!(τ_sub::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, τ_2::FT, f_sife::FT, τ_sub_sif::SubArray, vec_b::SubArray, vec_f::SubArray, ϕ::SubArray, N::Int) where {FT}

Update the SIF conversion matrix of the n-1 layer, given
- `τ_sub` transmittance within a sublayer
- `τ_θ` transmittance of the incoming radiation across the leaf layer
- `ρ_1` reflectance of the first layer (1)
- `ρ_2` reflectance of the second layer (n-1)
- `τ_2` transmittance of the second layer (n-1)
- `f_sife` SIF excitation scaling factor (f_cab + f_car * ϕ_car)
- `τ_sub_sif` transmittance within a sublayer for SIF
- `vec_b` SIF vector backwards (SubArray of a Matrix)
- `vec_f` SIF vector forwards (SubArray of a Matrix)
- `ϕ` SIF emission PDF (SubArray of a Vector)
- `N` number of sublayers of each layer

"""
function layer_2_sif_vec!(τ_sub::FT, τ_θ::FT, ρ_1::FT, ρ_2::FT, τ_2::FT, f_sife::FT, τ_sub_sif::SubArray, vec_b::SubArray, vec_f::SubArray, ϕ::SubArray, N::Int) where {FT}
    # parameters required for the calculation
    α_sub = 1 - τ_sub;
    τ_all = τ_sub ^ N;

    # clear the SIF vector
    vec_b .= 0;
    vec_f .= 0;

    # 1. here we consider the n-1 layers as one single layer, and the SIF transmission within this effective layer is same as the computed τ_sub_2
    #    then we need to rescale the interface ρ and τ for the effective layer so that the computed layer level ρ and τ are same as computed
    ρ_i_12 = effective_ρ_12(ρ_2, τ_2, τ_all);
    ρ_i_21 = effective_ρ_21(ρ_2, τ_2, τ_all);
    τ_i_12 = 1 - ρ_i_12;

    # 2. besides the the direction transmission from the first layer, we also need to account for the radiation reflected by the upper layer (the denominator 1 - ρ_1 * ρ_2)
    #    after accounting for the reflection and transmission within the layer
    #        the total radiation that transmits through is τ_θ / (1 - ρ_1 * ρ_2) * τ_2
    #        the total radiation that reflects back is τ_θ / (1 - ρ_1 * ρ_2) * ρ_1
    #    combine with that of the first layer, the total radiation is conserved
    #        the radiation that goes into the 1st layer is τ_θ * ρ_2 / (1 - ρ_1 * ρ_2) * (1 - ρ_1)
    #        the radiation that goes into the n-1 layer is τ_θ / (1 - ρ_1 * ρ_2) * (1 - ρ_2)
    #        the sum is then τ_θ (conserved)
    #    now the radiation goes from up to down
    rad_i = τ_θ / (1 - ρ_1 * ρ_2) * τ_i_12 * f_sife;
    for i in 1:N
        vec_b .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 3. then the radiation goes from down to up
    rad_i *= ρ_i_21;
    for i in N:-1:1
        vec_b .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (i - FT(0.5));
        vec_f .+= rad_i / 2 * α_sub * ϕ .* τ_sub_sif .^ (N - i + FT(0.5));
        rad_i *= τ_sub;
    end;

    # 4. rescale it by 1 / (1 - ρ_21 * t_all * ρ_21 * t_all) for the SIFE
    denom = 1 - τ_all * ρ_i_21 * τ_all * ρ_i_21;
    vec_b ./= denom;
    vec_f ./= denom;

    return nothing
end;
