aci_fit(config::SPACConfig{FT},
        ps::LeafPhotosystem{FT},
        pst::C4Trait{FT},
        acm::AcMethodC4Vcmax,
        ajm::AjMethodC4JPSII,
        apm::ApMethodC4VpmaxPi,
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
        if "Vpmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            iguess = findfirst(params .== "Vpmax25");
            push!(x_inis, initial_guess[iguess]);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
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
        if "Vcmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            push!(x_inis, 50);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Vpmax25" in params
            push!(x_mins, 1);
            push!(x_maxs, 200);
            push!(x_inis, 50);
            push!(Δ_inis, 10);
            push!(Δ_tols, 0.1);
        end;
        if "Rd25" in params
            push!(x_mins, 0.1);
            push!(x_maxs, 10);
            push!(x_inis, 1);
            push!(Δ_inis, 1);
            push!(Δ_tols, 0.01);
        end;
    end;

    mthd = ReduceStepMethodND{FT}(x_mins = x_mins, x_maxs = x_maxs, x_inis = x_inis, Δ_inis = Δ_inis);
    stol = SolutionToleranceND{FT}(Δ_tols, 50);
    func(x) = -aci_rmse(config, ps, pst, air, df, params, x);
    sol = find_peak(func, mthd, stol);

    best_rmse = aci_rmse(config, ps, pst, air, df, params, sol);
    aci = aci_curve(config, ps, air, df);

    return sol, best_rmse, aci
);


aci_rmse(config::SPACConfig{FT},
         ps::LeafPhotosystem{FT},
         pst::C4Trait{FT},
         acm::AcMethodC4Vcmax,
         ajm::AjMethodC4JPSII,
         apm::ApMethodC4VpmaxPi,
         air::AirLayer{FT},
         df::DataFrame,
         params::Vector{String},
         xxx::Vector) where {FT} = (
    if "Vcmax25" in params
        iparam = findfirst(params .== "Vcmax25");
        pst.v_cmax25 = xxx[iparam];
    end;
    if "Vpmax25" in params
        iparam = findfirst(params .== "Vpmax25");
        pst.v_pmax25 = xxx[iparam];
    end;
    if "Rd25" in params
        iparam = findfirst(params .== "Rd25");
        pst.r_d25 = xxx[iparam];
    end;

    return rmse(aci_curve(config, ps, air, df), df.A_NET)
);
