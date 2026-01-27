"""

$(TYPEDEF)

Struct that contains the trait variables for C3 photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct GeneralC3Trait{FT}
    # Related to OCS uptake
    "Multiplier to derive internal conductance for OCS `[mol μmol⁻¹]`"
    K_OCS::FT = 1400 * 1e-6

    # Prognostic variables
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT = 350 / 300
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    j_max25::FT = 83.5
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
end;


"""

$(TYPEDEF)

Struct that contains the trait variables for C4 photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct GeneralC4Trait{FT}
    # Related to OCS uptake
    "Multiplier to derive internal conductance for OCS `[mol μmol⁻¹]`"
    K_OCS::FT = 8862 * 1e-6

    # Prognostic variables
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_pmax25::FT = 50
end;


"""

$(TYPEDEF)

Struct that contains the state variables for C3 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3State{FT}
    # General model information
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT = 4
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT = 8

    # Prognostic variables (for VJP model)
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    k_npq_sus::FT = 0
end;


"""

$(TYPEDEF)

Struct that contains the state variables for C4 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C4State{FT}
    # Prognostic variables
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    k_npq_sus::FT = 0
end;


"""

$(TYPEDEF)

Struct that contains the auxiliary variables for leaf photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafPhotosystemAuxil{FT}
    # photosynthetic rates
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_c::Vector{FT}
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_g::Vector{FT}
    "Intermediate photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_i::Vector{FT}
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_j::Vector{FT}
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_n::Vector{FT}
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_p::Vector{FT}

    # electron transport rates
    "Electron to CO₂ coefficient"
    e2c::Vector{FT}
    "Fraction of absorbed light used by PSII ETR"
    f_psii::FT = 0.5
    "Electron transport `[μmol m⁻² s⁻¹]`"
    j::Vector{FT}
    "Maximal electron transport rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    j_max::FT = 0
    "Potential Electron Transport Rate `[μmol m⁻² s⁻¹]`"
    j_pot::Vector{FT}
    "PSI electron transport rate after colimitation"
    j_psi::Vector{FT}

    # photosynthesis rate coefficients
    "RubisCO coefficient Kc `[Pa]`"
    k_c::FT = 0
    "Michaelis-Menten's coefficient `[Pa]`"
    k_m::FT = 0
    "RubisCO coefficient Ko `[Pa]`"
    k_o::FT = 0
    "PEP coefficient Kpep `[Pa]`"
    k_pep::FT = 0
    "PEP coefficient Kpep fro CLM (different algorithm) `[Pa]`"
    k_pep_clm::FT = 0
    "Maximal turnover rate of Cytochrome b₆f `[e⁻ s⁻¹]`"
    k_q::FT = 0
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    γ_star::FT = 0

    # respiration and carboxylation
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    r_d::FT = 0
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_cmax::FT = 0
    "Maximal PEP carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_pmax::FT = 0
    "Maximal Cytochrome b₆f activity `[μmol e⁻ m⁻² s⁻¹]`"
    v_qmax::FT = 0

    # C3 Cytochrome model variables
    "ratio between J_P700 and J_P680"
    η::Vector{FT}
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_c::FT = 0
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_l::FT = 0

    # yield variables
    "Heat dissipation yield"
    ϕ_d::Vector{FT}
    "Fluorescence yield"
    ϕ_f::Vector{FT}
    "Non-photochemical yield"
    ϕ_n::Vector{FT}
    "Photochemical yield"
    ϕ_p::Vector{FT}

    # fluorescence yeolds of two photosystems
    "Fluorescence yield of PSI"
    ϕ_f1::Vector{FT}
    "Fluorescence yield of PSII"
    ϕ_f2::Vector{FT}

    # fluorescence variables
    "Dark adapted yield (`Kp=0`)"
    f_m::FT = 0
    "Light adapted yield (`Kp=0`)"
    f_m′::Vector{FT}
    "Dark-adapted fluorescence yield (`Kp=max`)"
    f_o::FT = 0
    "Light-adapted fluorescence yield in the dark (`Kp=max`)"
    f_o′::Vector{FT}
    "Non-Photochemical quenching "
    npq::Vector{FT}
    "Energy quenching"
    q_e::Vector{FT}
    "Photochemical quenching"
    q_p::Vector{FT}

    # fluorescence rate coefficients
    "Rate constant for thermal dissipation"
    k_d::FT = 0
    "Reversible NPQ rate constant (initially zero)"
    k_n::Vector{FT}
    "Rate constant for photochemistry"
    k_p::Vector{FT}
    "Maximal PS I photochemical yield"
    ϕ_psi_max::FT = 0
    "max PSII yield (_k_npq_rev = 0, all RC open)"
    ϕ_psii_max::FT = 0
end;

LeafPhotosystemAuxil(config::SPACConfig{FT}) where {FT} = (
    cache_dim_ppar = isnothing(config.DIMENSIONS.DIM_PPAR_BINS) ? config.DIMENSIONS.DIM_INCL * config.DIMENSIONS.DIM_AZI : config.DIMENSIONS.DIM_PPAR_BINS;

    return LeafPhotosystemAuxil{FT}(cache_dim_ppar+1)
);

LeafPhotosystemAuxil{FT}(dim::Int) where {FT} = (
    return LeafPhotosystemAuxil{FT}(
                a_c   = zeros(FT, dim),
                a_g   = zeros(FT, dim),
                a_i   = zeros(FT, dim),
                a_j   = zeros(FT, dim),
                a_n   = zeros(FT, dim),
                a_p   = zeros(FT, dim),
                e2c   = zeros(FT, dim),
                j     = zeros(FT, dim),
                j_pot = zeros(FT, dim),
                j_psi = zeros(FT, dim),
                η     = zeros(FT, dim),
                ϕ_d   = zeros(FT, dim),
                ϕ_f   = zeros(FT, dim),
                ϕ_n   = zeros(FT, dim),
                ϕ_p   = zeros(FT, dim),
                ϕ_f1  = zeros(FT, dim),
                ϕ_f2  = zeros(FT, dim),
                f_m′  = zeros(FT, dim),
                f_o′  = zeros(FT, dim),
                npq   = zeros(FT, dim),
                q_e   = zeros(FT, dim),
                q_p   = zeros(FT, dim),
                k_n   = zeros(FT, dim),
                k_p   = zeros(FT, dim)
    )
);


"""

$(TYPEDEF)

Struct that contains the fields for C3 photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafPhotosystem{FT}
    "Trait variables"
    trait::Union{GeneralC3Trait{FT}, GeneralC4Trait{FT}}
    "State variables"
    state::Union{C3State{FT}, C4State{FT}}
    "Auxilary variables"
    auxil::LeafPhotosystemAuxil{FT}
end;

LeafPhotosystem(config::SPACConfig{FT}, c3c4::String = "C3") where {FT} = (
    @assert c3c4 in ["C3", "C4"] "The model string should be either C3 or C4!";

    return if c3c4 == "C3"
        LeafPhotosystem{FT}(trait = GeneralC3Trait{FT}(), state = C3State{FT}(), auxil = LeafPhotosystemAuxil(config));
    else
        ps = LeafPhotosystem{FT}(trait = GeneralC4Trait{FT}(), state = C4State{FT}(), auxil = LeafPhotosystemAuxil(config));
        ps.auxil.f_psii = 0.41;
        ps
    end;
);
