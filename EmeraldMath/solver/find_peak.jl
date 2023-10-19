#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: move function to Emerald
#
#######################################################################################################################################################################################################
"""

    next_xy!(f::Function, xy::Matrix{FT}, history::Vector{Vector{FT}}, stepping::Bool) where {FT}

Determine the next points to simulate, given
- `f` Function to find peak
- `xy` Matrix of x (1st column) and y (2nd column)
- `history` A vector to save simulations
- `stepping` Optional. If true, save the optimization steps to the history field in method struct

"""
function next_xy!(f::Function, xy::Matrix{FT}, history::Vector{Vector{FT}}, stepping::Bool) where {FT}
    _x1,_x2,_x3,_y1,_y2,_y3 = xy;
    # if the curve if flat, do nothing
    if _y1==_y2==_y3
        return false

    # cases to update only the first part
    elseif (_y1==_y2>_y3) || (_y1>_y2)
        _xn = (_x1+_x2) / 2;
        _yn = f(_xn);
        if stepping
            push!(history, [_xn, _yn]);
        end;
        xy[3,1] = _x2;
        xy[3,2] = _y2;
        xy[2,1] = _xn;
        xy[2,2] = _yn;

        return true

    # cases to update only the second part
    elseif _y1<=_y2<_y3
        _xn = (_x2+_x3) / 2;
        _yn = f(_xn);
        if stepping
            push!(history, [_xn, _yn]);
        end;
        xy[1,1] = _x2;
        xy[1,2] = _y2;
        xy[2,1] = _xn;
        xy[2,2] = _yn;

        return true

    # cases to update both parts, _y1<_y2>=_y3
    else
        _xn = (_x1+_x2) / 2;
        _yn = f(_xn);
        if stepping
            push!(history, [_xn, _yn]);
        end;
        if _yn >= _y2
            xy[3,1] = _x2;
            xy[3,2] = _y2;
            xy[2,1] = _xn;
            xy[2,2] = _yn;
        else
            xy[1,1] = _xn;
            xy[1,2] = _yn;

            _xn = (_x2+_x3) / 2;
            _yn = f(_xn);
            if stepping
                push!(history, [_xn, _yn]);
            end;
            if _yn > _y2
                xy[1,1] = _x2;
                xy[1,2] = _y2;
                xy[2,1] = _xn;
                xy[2,2] = _yn;
            else
                xy[3,1] = _xn;
                xy[3,2] = _yn;
            end;
        end;

        return true
    end;
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-10: move function to Emerald
#
#######################################################################################################################################################################################################
"""

    find_peak(f::Function, ms::BisectionMethod{FT}, tol::Union{ResidualTolerance{FT}, SolutionTolerance{FT}}; stepping::Bool = false) where {FT}
    find_peak(f::Function, ms::NelderMeadMethod{FT}, tol::ResidualTolerance{FT}; stepping::Bool = false) where {FT}
    find_peak(f::Function, ms::NelderMeadMethod{FT}, tol::SolutionToleranceND{FT}; stepping::Bool = false) where {FT}
    find_peak(f::Function, ms::ReduceStepMethod{FT}, tol::SolutionTolerance{FT}; stepping::Bool = false) where {FT}
    find_peak(f::Function, ms::ReduceStepMethodND{FT}, tol::SolutionToleranceND{FT}; stepping::Bool = false) where {FT}

Find the solution that the y is maximum, given
- `f` A function to solve
- `ms` [`BisectionMethod`](@ref) type method struct
- `tol` [`SolutionTolerance`](@ref) type tolerance struct
- `stepping` Optional. If true, save the optimization steps to the history field in method struct.

"""
function find_peak end;

find_peak(f::Function, ms::BisectionMethod{FT}, tol::Union{ResidualTolerance{FT}, SolutionTolerance{FT}}; stepping::Bool = false) where {FT} = (
    # create matrix to store data
    (; history, x_max, x_min, xy) = ms;
    xy[1,1] = x_min;
    xy[2,1] = (x_min + x_max) / 2;
    xy[3,1] = x_max;
    xy[1,2] = f(xy[1,1]);
    xy[2,2] = f(xy[2,1]);
    xy[3,2] = f(xy[3,1]);
    if stepping
        push!(history, [xy[1,1], xy[1,2]]);
        push!(history, [xy[2,1], xy[2,2]]);
        push!(history, [xy[3,1], xy[3,2]]);
    end;

    # update the xy depending on the curve shape
    _count::Int = 0;
    while next_xy!(f, xy, history, stepping)
        _count += 1;
        if if_break(tol, xy[1,1], xy[3,1], xy[1,2]-xy[3,2], _count)
            break;
        end;
    end;
    (_maxy,_indy) = findmax( view(xy, :, 2) );

    return xy[_indy,1]
);

find_peak(f::Function, ms::NelderMeadMethod{FT}, tol::ResidualTolerance{FT}; stepping::Bool = false) where {FT} = (
    (; N, cen_x, con_x, exp_x, history, ref_x, simplex, x_inis) = ms;
    _nX = N;
    _nS = N + 1;

    # initialize the simplex
    for _irow in 1:_nS
        simplex[_irow] .= x_inis;
        if _irow < _nS
            simplex[_irow][_irow] = x_inis[_irow] / 2;
        end;
        simplex[_irow][end] = f(simplex[_irow]);
        if stepping
            push!(history, simplex[_irow]);
        end;
    end;

    # find the optimal simplex
    _count_all = 0;
    while true
        # 1. sort the simplex from high to low
        sort!(simplex, by=x->x[end], rev=true);
        if (simplex[1][end] - simplex[end][end]) < tol.tol
            break;
        end;

        # 2. calculate the centroid of the N best
        cen_x .= 0;
        for irow in 1:_nX
            cen_x .+= simplex[irow]
        end;
        cen_x ./= _nX;

        # 3. reflection of the worst
        ref_x  .= cen_x;
        ref_x .*= 2;
        ref_x .-= simplex[_nS];
        ref_x[end] = f(ref_x);
        if stepping
            push!(history, ref_x);
        end;

        if simplex[1][end] >= ref_x[end] > simplex[_nX][end]
            simplex[end] .= ref_x;

        # 4. expansion of the reflection
        elseif ref_x[end] > simplex[1][end]
            exp_x  .= ref_x;
            exp_x .*= 2;
            exp_x .-= cen_x;
            exp_x[end] = f(exp_x);
            if stepping
                push!(history, exp_x);
            end;

            if exp_x[end] > ref_x[end]
                simplex[end] .= exp_x;
            else
                simplex[end] .= ref_x;
            end;

        # 5. contraction of the worst
        else
            con_x  .= cen_x;
            con_x .+= simplex[_nS];
            con_x ./= 2;
            con_x[end] = f(con_x);
            if stepping
                push!(history, con_x);
            end;

            if con_x[end] > simplex[_nS][end]
                simplex[end] .= con_x;

            # 6. shrink
            else
                for _irow in 2:_nS
                    simplex[_irow] .+= simplex[1];
                    simplex[_irow] ./= 2;
                    simplex[_irow][end] = f(simplex[_irow]);
                    if stepping
                        push!(history, simplex[_irow]);
                    end;
                end;
            end;
        end;

        # 7. iteration ++
        _count_all += 1;

        # 8. determine whether to break
        _count_all > tol.n_limit ? break : nothing;
    end;

    return view(simplex[1], 1:_nX)
);

find_peak(f::Function, ms::NelderMeadMethod{FT}, tol::SolutionToleranceND{FT}; stepping::Bool = false) where {FT} = (
    (; cen_x, con_x, exp_x, history, ref_x, simplex, x_inis) = ms;
    _nX = ms.N;
    _nS = ms.N + 1;

    # initialize the simplex
    for _irow in 1:_nS
        simplex[_irow] .= x_inis
        if _irow < _nS
            simplex[_irow][_irow] = x_inis[_irow] / 2;
        end;
        simplex[_irow][end] = f(simplex[_irow])
        if stepping
            push!(history, simplex[_irow]);
        end;
    end;

    # find the optimal simplex
    _count_all = 0;
    judges = [false for i in 1:_nX];
    while true
        # 1. sort the simplex from high to low
        sort!(simplex, by=x->x[end], rev=true);
        for icol in 1:_nX
            judges[icol] = ( abs(simplex[1][icol] - simplex[_nS][icol]) <= tol.tol[icol]);
        end;
        all(judges) ? break : nothing;

        # 2. calculate the centroid of the N best
        cen_x .= FT(0)
        for irow in 1:_nX
            cen_x .+= simplex[irow]
        end;
        cen_x ./= _nX;

        # 3. reflection of the worst
        ref_x  .= cen_x;
        ref_x .*= 2;
        ref_x .-= simplex[_nS];
        ref_y   = f(ref_x);
        if stepping
            push!(history, ref_x);
        end;

        if simplex[1][end] >= ref_y > simplex[_nX][end]
            simplex[end]     .= ref_x;
            simplex[end][end] = ref_y;

        # 4. expansion of the reflection
        elseif ref_y > simplex[1][end]
            exp_x  .= ref_x;
            exp_x .*= 2;
            exp_x .-= cen_x;
            exp_y   = f(exp_x);
            if stepping
                push!(history, exp_x);
            end;

            if exp_y > ref_y
                simplex[end]     .= exp_x;
                simplex[end][end] = exp_y;
            else
                simplex[end]     .= ref_x;
                simplex[end][end] = ref_y;
            end;

        # 5. contraction of the worst
        else
            con_x  .= cen_x;
            con_x .+= simplex[_nS];
            con_x ./= 2;
            con_y   = f(con_x);
            if stepping
                push!(history, con_x);
            end;

            if con_y > simplex[_nS][end]
                simplex[end]     .= con_x;
                simplex[end][end] = con_y;

            # 6. shrink
            else
                for _irow in 2:_nS
                    simplex[_irow]    .+= simplex[1];
                    simplex[_irow][end] = f(simplex[_irow]);
                    if stepping
                        push!(history, simplex[_irow]);
                    end;
                end;
            end;
        end;

        # 7. iteration ++
        _count_all += 1;

        # 8. determine whether to break
        _count_all > tol.n_limit ? break : nothing;
    end;

    return view(simplex[1], 1:_nX)
);

find_peak(f::Function, ms::ReduceStepMethod{FT}, tol::SolutionTolerance{FT}; stepping::Bool = false) where {FT} = (
    # count iterations
    count_all = 0;

    # define the initial step
    (; x_ini, x_max, x_min, Δ_ini) = ms;

    # initialize the y
    tar_x = x_ini;
    tar_y = f(tar_x);
    new_x = x_min;
    new_y = f(new_x);
    new_y >= tar_y ? (tar_x=new_x; tar_y=new_y;) : nothing;
    new_x = x_max;
    new_y = f(new_x);
    new_y >  tar_y ? (tar_x=new_x; tar_y=new_y;) : nothing;

    # find the solution
    Δx::FT = Δ_ini;
    while true
        # 1. increase the x by Δx till tar_y is bigger
        count_inc = 0
        while true
            new_x = tar_x + Δx;
            new_x > x_max ? break : nothing;
            new_y = f(new_x);
            new_y > tar_y ? (tar_x=new_x; tar_y=new_y;) : break;
            count_inc += 1;
            count_all += 1;
        end;

        # 2. decrease the x by Δx till tar_y is bigger
        count_dec = 0
        while count_inc == 0
            new_x = tar_x - Δx;
            new_x < x_min ? break : nothing;
            new_y = f(new_x);
            new_y >= tar_y ? (tar_x=new_x; tar_y=new_y;) : break;
            count_dec += 1;
            count_all += 1;
        end;

        # 3. if no update, then 10% the Δx
        if count_inc + count_dec == 0
            Δx /= 10;
        end;

        # 4. if break
        Δx <= tol.tol ? break : nothing;
    end;

    return tar_x
);

find_peak(f::Function, ms::ReduceStepMethodND{FT}, tol::SolutionToleranceND{FT}; stepping::Bool = false) where {FT} = (
    # define the initial step
    (; x_maxs, x_mins, x_targ, x_temp, Δ_oper, Δjd) = ms;

    # initial the value
    tar_y = f(x_targ);

    # find the solution
    Nxs = length(x_targ);
    while true
        # total count of updates for all the x values
        count_all = 0;

        for ith in 1:Nxs
            # 1. increase the x by Δx till tar_y is bigger
            count_inc = 0
            while true
                x_temp[ith] = x_targ[ith] + Δ_oper[ith];
                x_temp[ith] > x_maxs[ith] ? break : nothing;
                new_y = f(x_temp);
                new_y >  tar_y ? (x_targ .= x_temp; tar_y=new_y) : break;
                count_inc += 1;
                count_all += 1;
            end;
            x_temp .= x_targ;

            # 2. decrease the x by Δx till tar_y is bigger
            while count_inc == 0
                x_temp[ith] = x_targ[ith] - Δ_oper[ith];
                x_temp[ith] < x_mins[ith] ? break : nothing;
                new_y = f(x_temp);
                new_y >= tar_y ? (x_targ .= x_temp; tar_y=new_y) : break;
                count_all += 1;
            end;
            x_temp .= x_targ;
        end;

        # 3. if no update, then 10% the Δx
        if count_all == 0
            for ith in 1:Nxs
                Δ_oper[ith] > tol.tol[ith] ? (Δ_oper[ith] /= 10) : (Δjd[ith] = true);
            end;
        end;

        # 4. judge whether to break
        all( Δjd ) ? break : nothing;
    end;

    return x_targ
);
