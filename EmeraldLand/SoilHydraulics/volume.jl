

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#    2023-Jun-30: move function out of soil_budget!
#     2023-Sep-07: add integrators for soil water budget
#
#######################################################################################################################################################################################################
"""

    surface_runoff!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}

Compute surface runoff, given
- `config` spac configuration
- `spac` spac model

"""
function surface_runoff!(config::SPACConfiguration{FT}, spac::MultiLayerSPAC{FT}) where {FT}
    (; SOIL_BULK, SOILS) = spac;

    # compute surface runoff
    if SOILS[1].state.θ > SOILS[1].state.vc.Θ_SAT
        # compute top soil temperature and top soil energy out due to runoff
        _cp_gas = (SOILS[1].state.ns[3] * CP_V_MOL(FT) + (SOILS[1].state.ns[1] + SOILS[1].state.ns[2] + SOILS[1].state.ns[4] + SOILS[1].state.ns[5]) * CP_D_MOL(FT)) / SOILS[1].auxil.δz;
        _cp = SOILS[1].state.cp * SOILS[1].state.ρ + SOILS[1].state.θ * ρ_H₂O(FT) * CP_L(FT) + _cp_gas;
        _t  = SOILS[1].state.Σe / _cp;
        _runoff = (SOILS[1].state.θ - SOILS[1].state.vc.Θ_SAT) * SOILS[1].auxil.δz * ρ_H₂O(FT) / M_H₂O(FT);

        SOILS[1].state.θ = SOILS[1].state.vc.Θ_SAT;
        SOILS[1].state.Σe -= _runoff / SOILS[1].auxil.δz * CP_L_MOL(FT) * _t;
        SOIL_BULK.auxil.runoff += _runoff;
    end;

    return nothing
end
