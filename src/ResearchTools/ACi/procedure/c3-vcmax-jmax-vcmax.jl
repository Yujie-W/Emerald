aci_fit(config::SPACConfig{FT},
        ps::LeafPhotosystem{FT},
        pst::GeneralC3Trait{FT},
        acm::AcMethodC3VcmaxPi,
        ajm::AjMethodC3JmaxPi,
        apm::ApMethodC3Vcmax,
        air::AirLayer{FT},
        df::DataFrame,
        params::Vector{String},
        initial_guess::Union{Nothing, Vector}) where {FT} = (
    # if initial_guess is not provided, derive the ranges from the data, else use the provided initial_guess (limit is twice of the initial_guess)
    x_mins = FT[];
    x_maxs = FT[];
    x_inis = FT[];
    Δ_inis = FT[];
    Δ_tols = FT[];
    if !isnothing(initial_guess)
        @assert length(initial_guess) == length(params);
        if "Vcmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            iguess = findfirst(params .== "Vcmax25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Jmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 400);
            iguess = findfirst(params .== "Jmax25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Γstar25" in params
            push!(x_mins, 1);
            push!(x_maxs, 10);
            iguess = findfirst(params .== "Γstar25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
        if "Rd25" in params
            push!(x_mins, 0.1);
            push!(x_maxs, 10);
            iguess = findfirst(params .== "Rd25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
    else
        # loop through the data once to get the respiration and Γ_star limits
        rd_lim_min = 0.1;
        γ_lim_max = 40;
        for dfr in eachrow(df)
            # estimate the minimum respiration rate limit
            rd_lim_min = nanmax([rd_lim_min, -dfr.A_NET / temperature_correction(config.METHODS.TD_R_C3, dfr.T_LEAF)]);
            # estimate the maximum Γ_star limit
            γ_lim_max = nanmin([γ_lim_max, dfr.P_I / temperature_correction(config.METHODS.TD_Γ, dfr.T_LEAF)]);
        end;
        # loop through the data once again to guess the Vcmax and Jmax
        vcmax_guess = 5;
        jmax_guess = 10;
        ps.trait.r_d25 = rd_lim_min * 1.2;
        config.METHODS.TD_Γ.VAL_REF = γ_lim_max * 0.8;
        for dfr in eachrow(df)
            photosystem_temperature_dependence!(config, ps, air, dfr.T_LEAF);
            vcmax = (dfr.A_NET + ps.auxil.r_d) * (dfr.P_I + ps.auxil.k_m) / (dfr.P_I - ps.auxil.γ_star) / temperature_correction(config.METHODS.TD_VCMAX_C3, dfr.T_LEAF);
            vcmax_guess = nanmax([vcmax_guess, vcmax]);
            vcmax_guess = nanmin([vcmax_guess, 100]);
            jmax = (dfr.A_NET + ps.auxil.r_d) * (4*dfr.P_I + 8*ps.auxil.γ_star) / (dfr.P_I - ps.auxil.γ_star) * 1.2 / temperature_correction(config.METHODS.TD_JMAX, dfr.T_LEAF);
            jmax_guess = nanmax([jmax_guess, jmax]);
            jmax_guess = nanmin([jmax_guess, 200]);
        end;
        # set the initial guess
        if "Vcmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, vcmax_guess * 2.0);
            push!(x_inis, vcmax_guess * 1.0);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Jmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, jmax_guess * 2.0);
            push!(x_inis, jmax_guess * 1.0);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Γstar25" in params
            push!(x_mins, 1);
            push!(x_maxs, γ_lim_max * 1.0);
            push!(x_inis, γ_lim_max * 0.8);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
        if "Rd25" in params
            push!(x_mins, rd_lim_min * 1.0);
            push!(x_maxs, 5);
            push!(x_inis, rd_lim_min * 1.2);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
    end;

    mthd = ReduceStepMethodND{FT}(x_mins = x_mins, x_maxs = x_maxs, x_inis = x_inis, Δ_inis = Δ_inis);
    stol = SolutionToleranceND{FT}(Δ_tols, 50);
    # func(x) = (rme = aci_rmse(config, ps, pst, air, df, x); @info "C3VJP model" x rme ; -rme);
    func(x) = -aci_rmse(config, ps, pst, air, df, params, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, params, sol);
    aci = aci_curve(config, ps, air, df);

    return sol, best_rmse, aci
);


aci_rmse(config::SPACConfig{FT},
         ps::LeafPhotosystem{FT},
         pst::GeneralC3Trait{FT},
         acm::AcMethodC3VcmaxPi,
         ajm::AjMethodC3JmaxPi,
         apm::ApMethodC3Vcmax,
         air::AirLayer{FT},
         df::DataFrame,
         params::Vector{String},
         xxx::Vector) where {FT} = (
    if "Vcmax25" in params
        iparam = findfirst(params .== "Vcmax25");
        pst.v_cmax25 = xxx[iparam];
    end;
    if "Jmax25" in params
        iparam = findfirst(params .== "Jmax25");
        pst.j_max25 = xxx[iparam];
    end;
    if "Γstar25" in params
        iparam = findfirst(params .== "Γstar25");
        config.METHODS.TD_Γ.VAL_REF = xxx[iparam];
    end;
    if "Rd25" in params
        iparam = findfirst(params .== "Rd25");
        pst.r_d25 = xxx[iparam];
    end;

    return rmse(aci_curve(config, ps, air, df), df.A_NET)
);
