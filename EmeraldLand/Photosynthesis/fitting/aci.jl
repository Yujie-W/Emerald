#######################################################################################################################################################################################################
#
# Changes to this functions
# General
#     2024-Jul-19: add functions to obtain net photosynthetic rate
#
#######################################################################################################################################################################################################
"""

    aci_an(ps::LeafPhotosystem, air::AirLayer, p_i::Number, ppar::Number, t::Number)

Compute the net photosynthetic rate, given
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `p_i` Internal CO₂ partial pressure in `Pa`
- `ppar` Photosynthetic photon flux density in `µmol m⁻² s⁻¹`
- `t` Leaf temperature in `K`

"""
function aci_an(ps::LeafPhotosystem, air::AirLayer, p_i::Number, ppar::Number, t::Number)
    photosynthesis!(ps, air, p_i, ppar, t);

    return ps.auxil.a_n
end;


#######################################################################################################################################################################################################
#
# Changes to this functions
# General
#     2024-Jul-19: add functions to obtain net photosynthetic rates of given data
#
#######################################################################################################################################################################################################
"""

    aci_curve(ps::LeafPhotosystem, air::AirLayer, pis::Vector, ppars::Vector, ts::Vector)
    aci_curve(ps::LeafPhotosystem, air::AirLayer, df::DataFrame)

Compute the net photosynthetic rates, given
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct
- `pis` Internal CO₂ partial pressure in `Pa`
- `ppars` Photosynthetic photon flux density in `µmol m⁻² s⁻¹`
- `ts` Leaf temperature in `K`
- `df` DataFrame with columns `P_I`, `PPAR`, and `T_LEAF`

"""
function aci_curve end;

aci_curve(ps::LeafPhotosystem, air::AirLayer, pis::Vector, ppars::Vector, ts::Vector) = aci_an.((ps,), (air,), pis, ppars, ts);

aci_curve(ps::LeafPhotosystem, air::AirLayer, df::DataFrame) = aci_an.(ps, air, df.P_I, df.PPAR, df.T_LEAF);
