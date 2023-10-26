# This file contains files for the bulk state and auxiliary variables for the soil

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-05: add struct SoilBulkAuxil
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for soil bulk state variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilBulkState{FT}
    "Total area of the soil `[m²]`"
    area::FT = 500
    "Color class as in CLM"
    color::Int = 1
end;


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-05: add struct SoilBulkAuxil
#     2023-Oct-07: add field dndt
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for soil bulk auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilBulkAuxil{FT}
    # soil energy budget
    "Net diffuse radiation at top soil `[mW m⁻² nm⁻¹]`"
    e_net_dif::Vector{FT}
    "Net direct radiation at top soil `[mW m⁻² nm⁻¹]`"
    e_net_dir::Vector{FT}
    "Net longwave energy absorption `[W m⁻²]`"
    r_net_lw::FT = 0
    "Net shortwave energy absorption `[W m⁻²]`"
    r_net_sw::FT = 0
    "Weights of the four characteristic curves"
    weight::Vector{FT} = zeros(FT, 4)
    "Reflectance for longwave radiation"
    ρ_lw::FT = 0.06
    "Reflectance for shortwave radiation"
    ρ_sw::Vector{FT}

    # the effective rate among soil layers
    "Diffusion rate between layers per area `[mol m⁻² s⁻¹]`"
    dndt::Matrix{FT}
    "Soil hydraulic conductance between layers per area `[mol m⁻² s⁻¹ MPa⁻¹]`"
    k::Vector{FT}
    "Flux between layers per area `[mol m⁻² s⁻¹]`"
    q::Vector{FT}
    "Thermal flux between layers per area `[mol m⁻² s⁻¹]`"
    q_layers::Vector{FT}
    "Surface runoff due to heavy precipitation during the time step `[mol m⁻²]`"
    runoff::FT = 0
    "Soil temperature difference between layers `[MPa]`"
    δt::Vector{FT}
    "Soil metric potential difference between layers `[MPa]`"
    δψ::Vector{FT}
    "Soil thermal conductance between layers per area `[W m⁻² K⁻¹]`"
    λ_layers::Vector{FT}

    # cache variables
    "Last soil moisture used to compute albedo"
    _θ::FT = -1
end;

SoilBulkAuxil(config::SPACConfiguration{FT}, n_soil::Int) where {FT} = SoilBulkAuxil{FT}(
            e_net_dif = zeros(FT, length(config.SPECTRA.Λ)),
            e_net_dir = zeros(FT, length(config.SPECTRA.Λ)),
            ρ_sw      = zeros(FT, length(config.SPECTRA.Λ)),
            dndt      = zeros(FT, n_soil, 5),
            k         = zeros(FT, n_soil - 1),
            q         = zeros(FT, n_soil - 1),
            q_layers  = zeros(FT, n_soil - 1),
            δt        = zeros(FT, n_soil - 1),
            δψ        = zeros(FT, n_soil - 1),
            λ_layers  = zeros(FT, n_soil - 1)
);


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Oct-05: add struct SoilBulk
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for soil bulk state and auxiliary variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SoilBulk{FT}
    "State variables"
    state::SoilBulkState{FT} = SoilBulkState{FT}()
    "Auxiliary variables"
    auxil::SoilBulkAuxil{FT}
end;

SoilBulk(config::SPACConfiguration{FT}, n_soil::Int) where {FT} = SoilBulk{FT}(auxil = SoilBulkAuxil(config, n_soil));
