leaf_sif_matrices_chl!(config::SPACConfig{FT}, bio::LeafBio{FT}, cache::SPACCache{FT}, ::DualspectFluorescenceSpectra) where {FT} = (
    (; SPECTRA) = config.CONSTANTS;
    (; IΛ_SIF, IΛ_SIFE, Λ_SIF, Λ_SIFE, Φ_PS) = SPECTRA;

    # update the SIF emission vector per excitation wavelength
    ϕ      = bio.auxil._ϕ_sif;
    factor = cache.cache_sif_1;

    for i in eachindex(IΛ_SIFE)
        ii = IΛ_SIFE[i];

        # read the SIF emission spectrum
        ϕ .= view(Φ_PS, IΛ_SIF);

        # tune SIF emission PDF based on the SIF excitation wavelength
        expsife = exp(Λ_SIFE[ii] / 10);
        @. factor = 1 / (1 + exp(-Λ_SIF / 10) * expsife);
        ϕ .*= factor;

        # read in the values from the auxiliary variables
        vec_b = view(bio.auxil.mat_b, :, i);
        vec_f = view(bio.auxil.mat_f, :, i);
        vec_b_chl = view(bio.auxil.mat_b_chl, :, i);
        vec_f_chl = view(bio.auxil.mat_f_chl, :, i);

        # compute the total absorbed radiation
        vec_b_chl .= bio.auxil.α_leaf[ii] * bio.auxil.f_sife[ii] .* ϕ .* vec_b ./ (vec_b .+ vec_f);
        vec_f_chl .= bio.auxil.α_leaf[ii] * bio.auxil.f_sife[ii] .* ϕ .* vec_f ./ (vec_b .+ vec_f);
    end;

    # compute the mean and mean diff of mat_b_chl and mat_f_chl
    bio.auxil.matꜛ_chl .= (bio.auxil.mat_b_chl .+ bio.auxil.mat_f_chl) ./ 2;
    bio.auxil.matꜜ_chl .= (bio.auxil.mat_b_chl .- bio.auxil.mat_f_chl) ./ 2;

    return nothing
);
