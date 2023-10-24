#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-09: add function to compute the weighted average of ϕ_d, ϕ_f, ϕ_n, and ϕ_p
#
#######################################################################################################################################################################################################
"""

    ϕDFNP(spac::BulkSPAC{FT}) where {FT}

Return the weighted average of ϕ_d, ϕ_f, ϕ_n, and ϕ_p, given
- `spac` `BulkSPAC` type struct

"""
function ΦDFNP(spac::BulkSPAC{FT}) where {FT}
    canopy = spac.canopy;
    leaves = spac.plant.leaves;

    sum_ϕda::FT = 0;
    sum_ϕfa::FT = 0;
    sum_ϕna::FT = 0;
    sum_ϕpa::FT = 0;
    sum_a::FT = 0;
    N = length(leaves);
    for i in eachindex(leaves)
        j = N - i + 1;
        sum_ϕda += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit .* leaves[i].flux.auxil.ϕ_d_sunlit) +
                   (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded * leaves[i].flux.auxil.ϕ_d_shaded) * canopy.structure.state.δlai[j];
        sum_ϕfa += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit .* leaves[i].flux.auxil.ϕ_f_sunlit) +
                   (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded * leaves[i].flux.auxil.ϕ_f_shaded) * canopy.structure.state.δlai[j];
        sum_ϕna += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit .* leaves[i].flux.auxil.ϕ_n_sunlit) +
                   (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded * leaves[i].flux.auxil.ϕ_n_shaded) * canopy.structure.state.δlai[j];
        sum_ϕpa += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit .* leaves[i].flux.auxil.ϕ_p_sunlit) +
                   (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded * leaves[i].flux.auxil.ϕ_p_shaded) * canopy.structure.state.δlai[j];
        sum_a   += (canopy.sun_geometry.auxil.p_sunlit[j] * mean(leaves[i].flux.auxil.a_g_sunlit) +
                   (1 - canopy.sun_geometry.auxil.p_sunlit[j]) * leaves[i].flux.auxil.a_g_shaded) * canopy.structure.state.δlai[j];
    end;

    return (sum_ϕda, sum_ϕfa, sum_ϕna, sum_ϕpa) ./ sum_a
end;
