
#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-13: add state struct to save
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Struct for states of monospecies tree SPAC system

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct MultiLayerSPACState{FT}
    # state variables
    "Shaded leaf stomatal conductance for all layers"
    gs_shaded::Vector{FT}
    "Sunlit leaf stomatal conductance for all layers"
    gs_sunlit::Array{FT,3}
    "Temperature record for CLM T mean of 10 days (based on CLM setting)"
    t_clm::Vector{FT}

    # variables to save
    "Soil moisture tuning factor β"
    beta::FT = 0
    "Total chloroplast fluorescence"
    csif::FT = 0
    "Total ETR"
    etr::FT = 0
    "Gross primary productivity"
    gpp::FT = 0
    "MODIS EVI"
    modis_evi::FT = 0
    "MODIS NDVI"
    modis_ndvi::FT = 0
    "MODIS NIRv"
    modis_nirv::FT = 0
    "OCO SIF at 759 nm"
    oco_sif₇₅₉::FT = 0
    "OCO SIF at 770 nm"
    oco_sif₇₇₀::FT = 0
    "PAR"
    par::FT = 0
    "PPAR"
    ppar::FT = 0
    "Transpiration"
    transpiration::FT = 0
    "TROPOMI SIF at 683 nm"
    tropomi_sif₆₈₃::FT = 0
    "TROPOMI SIF at 740 nm"
    tropomi_sif₇₄₀::FT = 0
end;

MultiLayerSPACState{FT}(spac::BulkSPAC{FT}) where {FT} = (
    leaves = spac.plant.leaves;

    gs_sunlit = zeros(FT, size(leaves[1].g_H₂O_s_sunlit,1), size(leaves[1].g_H₂O_s_sunlit,2), length(leaves));
    for i in eachindex(leaves)
        gs_sunlit[:,:,i] .= leaves[i].g_H₂O_s_sunlit;
    end;

    return MultiLayerSPACState{FT}(
                gs_shaded = [leaf.g_H₂O_s_shaded for leaf in leaves],
                gs_sunlit = gs_sunlit,
                t_clm = deepcopy(spac.plant.memory.t_history),
    )
);
