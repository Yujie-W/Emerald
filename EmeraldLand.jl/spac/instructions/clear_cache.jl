# This function contains functions to clear the cache in each time stepper iteration

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-30: add clear_∂X∂t! function
#
#######################################################################################################################################################################################################
"""

    clear_∂X∂t!(spac::MultiLayerSPAC{FT}) where {FT}

Clear the cache of ∂X∂t in each component of the SPAC, given
- `spac` `MultiLayerSPAC` SPAC

"""
function clear_∂X∂t!(spac::MultiLayerSPAC{FT}) where {FT}
    (; SOIL, ROOTS, JUNCTION, TRUNK, BRANCHES, LEAVES) = spac;

    # clear the dXdt cache in soil
    for soil in SOIL.LAYERS
        soil.∂e∂t = 0;
        soil.∂n∂t .= 0;
        soil.∂θ∂t = 0;
    end;

    # clear the dXdt cache in roots
    for root in ROOTS
        root.NS.energy.auxil.∂e∂t = 0;
    end;

    # clear the dXdt cache in junction
    JUNCTION.auxil.∂e∂t = 0;
    JUNCTION.auxil.∂w∂t = 0;

    # clear the dXdt cache in trunk
    TRUNK.NS.energy.auxil.∂e∂t = 0;

    # clear the dXdt cache in branches
    for branch in BRANCHES
        branch.NS.energy.auxil.∂e∂t = 0;
    end;

    # clear the dXdt cache in leaves
    for leaf in LEAVES
        leaf.NS.energy.auxil.∂e∂t = 0;
    end;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-30: add clear_∫∂X∂t! function
#
#######################################################################################################################################################################################################
"""

    clear_∫∂X∂t!(spac::MultiLayerSPAC{FT}) where {FT}

Clear the cache of ∫∂X∂t in each component of the SPAC, given
- `spac` `MultiLayerSPAC` SPAC

"""
function clear_∫∂X∂t!(spac::MultiLayerSPAC{FT}) where {FT}
    (; LEAVES) = spac;

    for leaf in LEAVES
        leaf.∫∂w∂t_out = 0;
    end;

    return nothing
end;
