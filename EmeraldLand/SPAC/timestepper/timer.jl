# This file contains function to calculate the time step size

#######################################################################################################################################################################################################
#
# Changes to the function
# General
#     2022-Jun-15: add function to make sure the soil layers do not over saturate or drain
#     2022-Jun-18: add controller for soil and leaf temperatures
#     2022-Aug-18: add option θ_on to enable/disable soil water budget
#     2022-Aug-31: add controller for leaf stomatal conductance
#     2022-Sep-07: remove soil oversaturation controller, and add a Δθ <= 0.01 controller
#     2022-Oct-22: add option t_on to enable/disable soil and leaf energy budgets
#     2023-Mar-27: add controller for trunk and branch temperatures
#     2023-Sep-30: add time controller of junction water
#     2023-Oct-18: add time controller of junction temperature
#
#######################################################################################################################################################################################################
"""

        adjusted_time(spac::BulkSPAC{FT}, δt::FT) where {FT}

Return adjusted time that soil does not over saturate or drain, given
- `config` Configuration for `BulkSPAC`
- `spac` `BulkSPAC` SPAC
- `δt` Time step

"""
function adjusted_time(spac::BulkSPAC{FT}, δt::FT) where {FT}
    branches = spac.plant.branches;
    junction = spac.plant.junction;
    lai = spac.canopy.structure.state.lai;
    leaves = spac.plant.leaves;
    soils = spac.soils;
    trunk = spac.plant.trunk;

    # make sure each layer does not drain (allow for oversaturation), and θ change is less than 0.01
    δt_1 = δt;
    for soil in soils
        δt_1 = min(FT(0.01) / abs(soil.auxil.∂θ∂t), δt_1);
        if soil.auxil.∂θ∂t < 0
            δt_drain = (soil.state.vc.Θ_RES - soil.state.θ) / soil.auxil.∂θ∂t;
            δt_1 = min(δt_drain, δt_1);
        end;
    end;

    # make sure temperatures do not change more than 1 K per time step
    δt_2 = δt_1;
    for soil in soils
        ∂T∂t = soil.auxil.∂e∂t / soil.auxil.cp;
        δt_2 = min(1 / abs(∂T∂t), δt_2);
    end;

    ∂T∂t = junction.auxil.∂e∂t / junction.auxil.cp;
    δt_2 = min(1 / abs(∂T∂t), δt_2);

    ∂T∂t = trunk.energy.auxil.∂e∂t / trunk.energy.auxil.cp;
    δt_2 = min(1 / abs(∂T∂t), δt_2);

    for stem in branches
        ∂T∂t = stem.energy.auxil.∂e∂t / stem.energy.auxil.cp;
        δt_2 = min(1 / abs(∂T∂t), δt_2);
    end;

    if lai > 0
        for leaf in leaves
            ∂T∂t = leaf.energy.auxil.∂e∂t / leaf.energy.auxil.cp;
            δt_2 = min(1 / abs(∂T∂t), δt_2);
        end;
    end;

    # make sure the junction water does not change more than 10 mol per time step
    δt_3 = δt_2;
    if junction.auxil.∂w∂t != 0
        δt_3 = min(10 / abs(junction.auxil.∂w∂t), δt_3);
    end;

    # make sure leaf stomatal conductances do not change more than 0.01 mol m⁻² s⁻¹
    δt_4 = δt_3;
    if lai > 0
        for leaf in leaves
            for ∂g∂t in leaf.flux.auxil.∂g∂t_sunlit
                δt_4 = min(FT(0.01) / abs(∂g∂t), δt_4);
            end;
            δt_4 = min(FT(0.01) / abs(leaf.flux.auxil.∂g∂t_shaded), δt_4);
        end;
    end;

    # make sure adjusted time is not nan
    if isnan(δt_4)
        @error "NaN detected in adjusted_time" δt δt_1 δt_2 δt_3 δt_4;
    end;

    return δt_4
end;
