# This file contains functions to compute leaf photosynthesis of the entire plant

######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-29: add method for BulkSPAC
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2023-Mar-11: only compute respiration rate if solar zenith angle >= 89
#     2023-Mar-11: do nothing if LAI == 0
#     2024-Jul-25: save average a_g and a_n (to use later)
#     2024-Jul-30: compute OCS fluxes along with photosynthesis
#
#######################################################################################################################################################################################################
"""

    plant_photosynthesis!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Updates leaf photosynthetic rates for SPAC, given
- `config` `SPACConfiguration` type SPAC configuration
- `spac` `BulkSPAC` type SPAC

"""
function plant_photosynthesis! end;

plant_photosynthesis!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT} = plant_photosynthesis!(config, spac, spac.plant.leaves[1]);

plant_photosynthesis!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}, ::CanopyLayer{FT}) where {FT} = (
    # if there is no leaf, do nothing
    if spac.canopy.structure.trait.lai <= 0
        return nothing
    end;

    airs = spac.airs;
    canopy = spac.canopy;
    leaves = spac.plant.leaves;
    lindex = spac.plant.leaves_index;
    n_layer = length(leaves);

    rd_only = spac.canopy.sun_geometry.state.sza > 89;
    for ilf in eachindex(leaves)
        irt = n_layer + 1 - ilf;
        leaf = leaves[ilf];
        air = airs[lindex[ilf]];
        leaf_photosynthesis!(config, spac.cache, leaf, air; rd_only = rd_only);

        # update the OCS flux
        leaf.flux.auxil.f_ocs .= leaf.flux.auxil.g_OCS .* air.s_aux.ps[6] ./ air.state.p_air .* FT(1e6);

        # update the average rates
        leaf.flux.auxil.a_g_mean   = leaf.flux.auxil.a_g'   * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
        leaf.flux.auxil.a_n_mean   = leaf.flux.auxil.a_n'   * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
        leaf.flux.auxil.f_ocs_mean = leaf.flux.auxil.f_ocs' * view(canopy.sun_geometry.auxil.ppar_fraction,:,irt);
    end;

    return nothing
);


######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2024-Aug-29: add function to update the carbon pool budget of the plant
#     2024-Sep-03: do not update the carbon pool if the plant is dead
#
#######################################################################################################################################################################################################
"""

    plant_carbon_budget!(spac::BulkSPAC{FT}, δt::FT) where {FT}

Update the carbon pool budget of the plant, given
- `spac` `BulkSPAC` SPAC
- `δt` Time step

"""
function plant_carbon_budget!(spac::BulkSPAC{FT}, δt::FT) where {FT}
    plant = spac.plant;

    # if the plant is dead, do nothing
    if plant.pool.c_pool <= 0
        return nothing
    end;

    # TODO update the carbon pool of the plant to account for root and stem respiration
    # TODO chemical energy change
    for r in plant.roots
        c_mol = r.xylem.trait.area * r.xylem.trait.l * r.xylem.trait.ρ * 1000 / 30; # mol C
        resp = temperature_corrected_value(r.xylem.trait.r_wood, r.energy.s_aux.t); # μmol mol⁻¹ s⁻¹
        f = resp * FT(1e-6) * c_mol * δt;
        plant.pool.c_pool -= f;
    end;

    # trunk respiration
    c_mol = plant.trunk.xylem.trait.area * plant.trunk.xylem.trait.l * plant.trunk.xylem.trait.ρ * 1000 / 30; # mol C
    resp = temperature_corrected_value(plant.trunk.xylem.trait.r_wood, plant.trunk.energy.s_aux.t); # μmol mol⁻¹ s⁻¹
    f = resp * FT(1e-6) * c_mol * δt;
    plant.pool.c_pool -= f;

    # branch respiration
    for s in plant.branches
        c_mol = s.xylem.trait.area * s.xylem.trait.l * s.xylem.trait.ρ * 1000 / 30; # mol C
        resp = temperature_corrected_value(s.xylem.trait.r_wood, s.energy.s_aux.t); # μmol mol⁻¹ s⁻¹
        f = resp * FT(1e-6) * c_mol * δt;
        plant.pool.c_pool -= f;
    end;

    # do nothing if LAI == 0; otherwise update the carbon budget of each leaf (or canopy layer)
    if spac.canopy.structure.trait.lai <= 0
        return nothing
    end;

    for l in plant.leaves
        f = l.flux.auxil.a_n_mean * FT(1e-6) * l.xylem.trait.area * δt;
        l.flux.auxil.∫∂c∂t_in += f;
        plant.pool.c_pool += f;
    end;

    return nothing
end;
