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
    "Isotropic radiation"
    ISO_RAD::Vector{FT}

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
    Σkx::Vector{Matrix{FT}}
    "Downward transmittance for isotropic radiation"
    τꜜ::Vector{FT}
    "Upward transmittance for isotropic radiation"
    τꜛ::Vector{FT}

    # cache variables
    "Cache variables with length of wavelength"
    cache_wl_1::Vector{FT}
    "Cache variables with length of diffuse angles"
    cache_dif_1::Vector{FT}
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
                ISO_RAD     = sind.(collect(FT, 0.5:1:89.5)) ./ sum(sind.(collect(FT, 0.5:1:89.5))),
                # auxiliary variables
                ρ_12        = zeros(length(df.WL)),
                ρ_21        = zeros(length(df.WL)),
                τ_12        = zeros(length(df.WL)),
                τ_21        = zeros(length(df.WL)),
                Σkx         = [ zeros(FT, length(df.WL), 90) for i in 1:layers ],
                τꜜ          = zeros(FT, length(df.WL)),
                τꜛ          = zeros(FT, length(df.WL)),
                # cache
                cache_wl_1  = zeros(FT, length(df.WL)),
                cache_dif_1 = zeros(FT, 90)
    );
);


# step 1: compute the reflectance and transmittance at the air-water interface
function interface_ρ_τ!(model::OceanWaterLayers{FT}) where {FT}
    @. model.ρ_12 = interface_isotropic_τ(FT(1), model.NR, FT(90));
    @. model.ρ_21 = interface_isotropic_τ(model.NR, FT(1), FT(90));
    @. model.τ_12 = 1 - model.ρ_12;
    @. model.τ_21 = 1 - model.ρ_21;

    return nothing
end;


# step 2: compute the extinction coefficient per wavelength, per angle, per layer
function layer_extinction_coefficient!(model::OceanWaterLayers{FT}) where {FT}
    for i in eachindex(model.Δz)
        h2o_v = model.Δz[i] * 100; # cm
        chl_v = model.chl[i] * model.Δz[i] * 1000 / 10000; # ug cm⁻²
        # loop through the angles from 0.5 to 89.5
        for i_dif in 1:90
            @. model.Σkx[i][:,i_dif] = model.K_H₂O * h2o_v + model.K_CAB * chl_v;
            @. model.Σkx[i][:,i_dif] /= cosd(i_dif - FT(0.5));
        end;
    end;

    return nothing
end;


# step 3: loop through the layers to compute the downward transmission for isotropic radiation
function downward_τ!(model::OceanWaterLayers{FT}) where {FT}
    # loop throught the wavelengths
    for i_wl in eachindex(model.NR)
        e_rad = model.cache_dif_1;
        @. e_rad = model.ISO_RAD;
        # loop through the layers from top to bottom
        for i_layer in eachindex(model.Δz)
            @. e_rad *= exp(-model.Σkx[i_layer][i_wl,:]);
        end;
        model.τꜜ[i_wl] = sum(e_rad);
    end;

    return nothing
end;


# step 4: loop through the layers to compute the upward transmission for isotropic radiation
function upward_τ!(model::OceanWaterLayers{FT}) where {FT}
    # loop throught the wavelengths
    for i_wl in eachindex(model.NR)
        e_rad = model.cache_dif_1;
        @. e_rad = model.ISO_RAD;
        # loop through the layers from bottom to top
        nlayer = length(model.Δz);
        for i_layer in eachindex(model.Δz)
            @. e_rad *= exp(-model.Σkx[nlayer+1-i_layer][i_wl,:]);
        end;
        model.τꜛ[i_wl] = sum(e_rad);
    end;

    return nothing
end;


end # module


#= testing the module


using Emerald
using Revise

# create the structure
model = Emerald.EmeraldOcean.OceanOptics.OceanWaterLayers{Float64}(100);
Emerald.EmeraldOcean.OceanOptics.interface_ρ_τ!(model);
Emerald.EmeraldOcean.OceanOptics.layer_extinction_coefficient!(model);
Emerald.EmeraldOcean.OceanOptics.downward_τ!(model);
Emerald.EmeraldOcean.OceanOptics.upward_τ!(model);
@show model.τꜜ;
@show model.τꜛ;


=#
