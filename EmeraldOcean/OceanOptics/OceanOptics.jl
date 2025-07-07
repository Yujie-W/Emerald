module OceanOptics

using LinearAlgebra: mul!

using ..EmeraldIO.Jld2: read_jld2
using ..EmeraldLand.Namespace: LAND_ARTIFACT, OLD_PHI_2021_1NM
using ..EmeraldLand.LeafOptics: interface_isotropic_τ


# struct for the model
Base.@kwdef mutable struct OceanWaterLayers{FT}
    # structure of ocean water layers
    "Layer thickness [m]"
    Δz::Vector{FT}
    "Chlorophyll content per layer [mg m⁻³]"
    chl::Vector{FT}
    "Bottom albedo per wavelength"
    ρ_bottom::Vector{FT}

    # chlorophyll fluorescence quantum yield
    "Chlorophyll fluorescence quantum yield per layer"
    ϕ_f::Vector{FT}

    # constants
    "Refractive index of the water per wavelength"
    NR::Vector{FT}
    "Absorption coefficient of water per wavelength"
    K_H₂O::Vector{FT}
    "Absorption coefficient of chlorophyll per wavelength"
    K_CAB::Vector{FT}
    "Wavelength"
    Λ::Vector{FT}
    "SIF spectrum"
    Φ_PS::Vector{FT}
    "Isotropic radiation"
    ISO_RAD::Vector{FT}
    "Isotropic SIF"
    ISO_SIF::Matrix{FT} = Φ_PS * ISO_RAD'

    # auxiliary variables to compute for convenience
    "Surface reflectance per wavelength at the air-water interface"
    ρ_12::Vector{FT}
    "Surface reflectance per wavelength at the water-air interface"
    ρ_21::Vector{FT}
    "Surface transmittance per wavelength at the air-water interface"
    τ_12::Vector{FT}
    "Surface transmittance per wavelength at the water-air interface"
    τ_21::Vector{FT}
    "Extinction coefficient per wavelength, per angle, per layer"
    Σkx::Array{FT,3}
    "Ratio of the absorption by chlorophyll per wavelength per layer"
    f_chl::Matrix{FT}
    "Downward transmittance for isotropic radiation"
    τꜜ::Vector{FT}
    "Upward transmittance for isotropic radiation"
    τꜛ::Vector{FT}
    "Absorption ratio per wavelength per layer for downward isotropic radiation"
    αꜜ::Matrix{FT}
    "Absorption ratio per wavelength per layer for upward isotropic radiation"
    αꜛ::Matrix{FT}

    # now the auxiliary variables to derive SIF
    "SIF emitted at the given layer per SIF wavelength, per angle, per layer (does not matter upward or downward, half each)"
    sif_chl_wle::Array{FT,3}
    "Downward SIF without any reflection per SIF wavelength, per angle"
    sifꜜ_wle::Matrix{FT}
    "Upward SIF without any reflection per SIF wavelength, per angle"
    sifꜛ_wle::Matrix{FT}
    "Downward SIF without any reflection per SIF wavelength"
    sifꜜ::Vector{FT}
    "Upward SIF without any reflection per SIF wavelength"
    sifꜛ::Vector{FT}
    "Escaped SIF from the surface"
    sif_esc::Vector{FT}

    # cache variables
    "Cache variables with length of wavelength"
    cache_wl_1::Vector{FT}
    "Cache variables with length of diffuse angles"
    cache_dif_1::Vector{FT}

    # outputs that can be observed by the user
    "Surface reflectance per wavelength"
    ρ_surface::Vector{FT}
    "SIF matrix"
    sif_mat::Matrix{FT}
end;


# constructor for the model
OceanWaterLayers{FT}(layers::Int; dataset::String = OLD_PHI_2021_1NM, jld2_file::String = LAND_ARTIFACT) where {FT} = (
    df = read_jld2(jld2_file, dataset);

    return OceanWaterLayers{FT}(
                # structure
                Δz          = ones(FT, layers) .* 0.1,
                chl         = zeros(FT,layers),
                ρ_bottom    = ones(length(df.WL)) .* 0.2,
                # SIF
                ϕ_f         = ones(FT, layers) .* 0.01,
                # constants
                NR          = df.NR,
                K_H₂O       = df.K_H₂O,
                K_CAB       = df.K_CAB,
                Λ           = df.WL,
                Φ_PS        = df.K_PS,
                ISO_RAD     = sind.(collect(FT, 0.5:1:89.5)) ./ sum(sind.(collect(FT, 0.5:1:89.5))),
                # auxiliary variables
                ρ_12        = zeros(length(df.WL)),
                ρ_21        = zeros(length(df.WL)),
                τ_12        = zeros(length(df.WL)),
                τ_21        = zeros(length(df.WL)),
                Σkx         = zeros(FT, length(df.WL), 90, layers),
                f_chl       = zeros(FT, length(df.WL), layers),
                τꜜ          = zeros(FT, length(df.WL)),
                τꜛ          = zeros(FT, length(df.WL)),
                αꜜ          = zeros(FT, length(df.WL), layers),
                αꜛ          = zeros(FT, length(df.WL), layers),
                # SIF auxiliary variables
                sif_chl_wle = zeros(FT, length(df.WL), 90, layers),
                sifꜜ_wle    = zeros(FT, length(df.WL), 90),
                sifꜛ_wle    = zeros(FT, length(df.WL), 90),
                sifꜜ        = zeros(FT, length(df.WL)),
                sifꜛ        = zeros(FT, length(df.WL)),
                sif_esc     = zeros(FT, length(df.WL)),
                # cache
                cache_wl_1  = zeros(FT, length(df.WL)),
                cache_dif_1 = zeros(FT, 90),
                # outputs
                ρ_surface   = zeros(length(df.WL)),
                sif_mat     = zeros(FT, length(df.WL), length(df.WL))
    );
);


# step 1: compute the reflectance and transmittance at the air-water interface
function interface_ρ_τ!(owl::OceanWaterLayers{FT}) where {FT}
    @. owl.τ_12 = interface_isotropic_τ(FT(1), owl.NR, FT(90));
    @. owl.τ_21 = interface_isotropic_τ(owl.NR, FT(1), FT(90));
    @. owl.ρ_12 = 1 - owl.τ_12;
    @. owl.ρ_21 = 1 - owl.τ_21;

    return nothing
end;


# step 2: compute the extinction coefficient per wavelength, per angle, per layer
function layer_extinction_coefficient!(owl::OceanWaterLayers{FT}) where {FT}
    for i in eachindex(owl.Δz)
        h2o_v = owl.Δz[i] * 100; # cm
        chl_v = owl.chl[i] * owl.Δz[i] * 1000 / 10000; # ug cm⁻²
        # loop through the angles from 0.5 to 89.5
        for i_dif in 1:90
            @. owl.Σkx[:,i_dif,i] = owl.K_H₂O * h2o_v + owl.K_CAB * chl_v;
            @. owl.Σkx[:,i_dif,i] /= cosd(i_dif - FT(0.5));
        end;
        @. owl.f_chl[:,i] = owl.K_CAB * chl_v / (owl.K_H₂O * h2o_v + owl.K_CAB * chl_v);
    end;

    return nothing
end;


# step 3: loop through the layers to compute the downward transmission for isotropic radiation
function downward_τ!(owl::OceanWaterLayers{FT}) where {FT}
    # loop throught the wavelengths
    for i_wl in eachindex(owl.NR)
        e_rad = owl.cache_dif_1;
        @. e_rad = owl.ISO_RAD;
        last_rad = 1;
        curr_rad = 1;
        # loop through the layers from top to bottom
        for i_layer in eachindex(owl.Δz)
            @. e_rad *= exp(-owl.Σkx[i_wl,:,i_layer]);
            curr_rad = sum(e_rad);
            owl.αꜜ[i_wl,i_layer] = (last_rad - curr_rad) * owl.f_chl[i_wl,i_layer];
            last_rad = curr_rad;
        end;
        owl.τꜜ[i_wl] = last_rad;
    end;

    return nothing
end;


# step 4: loop through the layers to compute the upward transmission for isotropic radiation
function upward_τ!(owl::OceanWaterLayers{FT}) where {FT}
    # loop throught the wavelengths
    for i_wl in eachindex(owl.NR)
        e_rad = owl.cache_dif_1;
        @. e_rad = owl.ISO_RAD;
        last_rad = 1;
        curr_rad = 1;
        # loop through the layers from bottom to top
        nlayer = length(owl.Δz);
        for i_layer in eachindex(owl.Δz)
            @. e_rad *= exp(-owl.Σkx[i_wl,:,nlayer+1-i_layer]);
            curr_rad = sum(e_rad);
            owl.αꜛ[i_wl,nlayer+1-i_layer] = (last_rad - curr_rad) * owl.f_chl[i_wl,nlayer+1-i_layer];
            last_rad = curr_rad;
        end;
        owl.τꜛ[i_wl] = last_rad;
    end;

    return nothing
end;


# step 5: compute the surface reflectance
function surface_ρ!(owl::OceanWaterLayers{FT}) where {FT}
    ρ_12 = owl.ρ_12;
    ρ_21 = owl.ρ_21;
    τ_12 = owl.τ_12;
    τ_21 = owl.τ_21;
    ρ_bm = owl.ρ_bottom;
    τꜜ   = owl.τꜜ;
    τꜛ   = owl.τꜛ;
    @. owl.ρ_surface = ρ_12 + τ_12 * τꜜ * ρ_bm * τꜛ * τ_21 / (1 - ρ_21 * τꜜ * ρ_bm * τꜛ);

    return nothing
end;


# step 6: compute the SIF emitted at each layer for a given wavelength of incoming radiation
function layer_sif!(owl::OceanWaterLayers{FT}, i_wle::Int) where {FT}
    # the total SIF emitted is the product of the quantum yield and the absorbed radiation
    for i in eachindex(owl.Δz)
        rad_e = owl.τ_12[i_wle] / (1 - owl.τꜜ[i_wle] * owl.ρ_bottom[i_wle] * owl.τꜛ[i_wle] * owl.ρ_21[i_wle]) * owl.αꜜ[i_wle,i] +
                owl.τ_12[i_wle] * owl.τꜜ[i_wle] * owl.ρ_bottom[i_wle] / (1 - owl.τꜜ[i_wle] * owl.ρ_bottom[i_wle] * owl.τꜛ[i_wle] * owl.ρ_21[i_wle]) * owl.αꜛ[i_wle,i];
        rad_f = rad_e * owl.ϕ_f[i];
        # half goes up, half goes down
        # account for reabsorption for half the length
        @. owl.sif_chl_wle[:,:,i] = owl.ISO_SIF * rad_f / 2;
        @. owl.sif_chl_wle[:,:,i] *= exp(-owl.Σkx[:,:,i] / 2);
    end;

    return nothing
end;


# step 7: compute the total SIF reaching the bottom
function interface_sif!(owl::OceanWaterLayers{FT}, i_wle::Int) where {FT}
    # loop through the layers from top to bottom to compute the SIF reaching the bottom
    @. owl.sifꜜ_wle = 0;
    for i_layer in eachindex(owl.Δz)
        @. owl.sifꜜ_wle *= exp(-owl.Σkx[:,:,i_layer]);
        @. owl.sifꜜ_wle += owl.sif_chl_wle[:,:,i_layer];
    end;

    # sum up the SIF at each angle
    @. owl.sifꜜ = 0;
    for i_dif in 1:90
        @. owl.sifꜜ += owl.sifꜜ_wle[:,i_dif];
    end;

    # loop through the layers from bottom to top to compute the SIF reaching the top
    @. owl.sifꜛ_wle = 0;
    layers = length(owl.Δz);
    for i_layer in layers:-1:1
        @. owl.sifꜛ_wle *= exp(-owl.Σkx[:,:,i_layer]);
        @. owl.sifꜛ_wle += owl.sif_chl_wle[:,:,i_layer];
    end;

    # sum up the SIF at each angle
    @. owl.sifꜛ = 0;
    for i_dif in 1:90
        @. owl.sifꜛ += owl.sifꜛ_wle[:,i_dif];
    end;

    return nothing
end;


# step 8: compute the SIF escaped from the surface
function escaped_sif!(owl::OceanWaterLayers{FT}, i_wle::Int) where {FT}
    @. owl.sif_esc = (owl.sifꜛ + owl.sifꜜ * owl.ρ_bottom * owl.τꜛ) * owl.τ_21 / (1 - owl.ρ_21 * owl.τꜜ * owl.ρ_bottom * owl.τꜛ);

    return nothing
end;


# step 9: compute the SIF matrix
#         total escaped SIF is owl.sif_mat * in_rad
#         the unit needs to be converted to photons prior to the computation, and then converted back to energy
function sif_matrix!(owl::OceanWaterLayers{FT}) where {FT}
    # loop through the incoming radiation wavelengths
    for i_wle in eachindex(owl.NR)
        if 300 <= owl.Λ[i_wle] <= 750
            @show owl.Λ[i_wle];
            layer_sif!(owl, i_wle);
            interface_sif!(owl, i_wle);
            escaped_sif!(owl, i_wle);
            @. owl.sif_mat[:,i_wle] = owl.sif_esc;
        end;
    end;

    return nothing
end;


end # module


#= testing the module

using DataFrames
using Emerald
using Revise

# create the structure
owl = Emerald.EmeraldOcean.OceanOptics.OceanWaterLayers{Float64}(100);
owl.Δz .= 1;
owl.chl[10] = 5;
@time begin
    Emerald.EmeraldOcean.OceanOptics.interface_ρ_τ!(owl);
    Emerald.EmeraldOcean.OceanOptics.layer_extinction_coefficient!(owl);
    Emerald.EmeraldOcean.OceanOptics.downward_τ!(owl);
    Emerald.EmeraldOcean.OceanOptics.upward_τ!(owl);
    Emerald.EmeraldOcean.OceanOptics.surface_ρ!(owl);
    Emerald.EmeraldOcean.OceanOptics.sif_matrix!(owl);
end;

# read the radiation
df = Emerald.EmeraldIO.Jld2.read_jld2(Emerald.EmeraldLand.Namespace.LAND_ARTIFACT, Emerald.EmeraldLand.Namespace.OLD_PHI_2021_1NM);
rad = df.E_DIR .+ df.E_DIFF;
rad_photon = Emerald.EmeraldPhysics.Optics.photon.(owl.Λ, rad);
sif_photon = owl.sif_mat * rad_photon;
sif = Emerald.EmeraldPhysics.Optics.energy.(owl.Λ, sif_photon);

df_owl = DataFrames.DataFrame(WL = owl.Λ, NR = owl.NR, REFL = owl.ρ_surface, SIF = sif, PHI = owl.Φ_PS);
Emerald.EmeraldIO.Text.save_csv!(df_owl, "sif.csv");


=#
