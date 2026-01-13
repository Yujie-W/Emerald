

import Emerald.LeafOptics as ELO
n₁ = 1.0;
n₂ = 1.4;

for θ_max in 1.0:90
    sum_sin_tau = 0.0;
    sum_sin = 0.0;
    for θ in 0.5:1:θ_max
        tau = ELO.interface_ρ_τ(n₁, n₂, θ)[1];
        sum_sin_tau += tau * sind(θ);
        sum_sin += sind(θ);
    end;
    @info "Debug" θ_max sum_sin_tau / sum_sin 1 - ELO.interface_isotropic_τ(n₁, n₂, θ_max);
end;



import Emerald.LeafOptics as ELO
import Emerald.Namespace as ENS
import Emerald.SPAC as ESPAC


config = ENS.SPACConfig(Float64);
spac = ENS.BulkSPAC(config);
ESPAC.initialize_spac!(config, spac);
lbio = spac.plant.leaves[end].bio;

ELO.leaf_spectra!(config, spac.plant.leaves[end].bio, spac.cache, 5.0);
lbio.auxil.ρ_layer_1[61]
lbio.auxil.τ_layer_1[61]

#=
function interface_isotropic_τ(n₁::FT, n₂::FT, θ₁::FT) where {FT}
    sum_sin_tau = 0;
    sum_sin = 0;
    for θ in FT(0.05):FT(0.1):θ₁
        tau = interface_ρ_τ(n₁, n₂, θ)[2];
        sum_sin_tau += tau * sind(θ);
        sum_sin += sind(θ);
    end;

    return sum_sin_tau / sum_sin
end;
=#
