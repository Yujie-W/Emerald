# This function contains functions to set up root flow profile

flow_in(root::Root{FT}) where {FT} = flow_in(root.xylem);

flow_out(root::Root{FT}) where {FT} = flow_out(root.xylem);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-23: add function root_flow_profile!
#
#######################################################################################################################################################################################################
"""

    root_flow_profile!(config::SPACConfiguration{FT}, root::Root{FT}, soil::SoilLayer{FT}, p_target::FT) where {FT}

Use solver to determine the root flow rate, given
- `config` `SPACConfiguration` type struct
- `root` `Root` type struct
- `soil` `SoilLayer` type struct
- `p_target` Target pressure at the end; of root xylem

"""
function root_flow_profile!(config::SPACConfiguration{FT}, root::Root{FT}, soil::SoilLayer{FT}, junction::JunctionCapacitor{FT}) where {FT}
    # if the root is not connected to the soil, set the flow to be the sum from the buffer system
    # else, use a solver to find the root flow rate
    if soil.auxil.ψ <= xylem_pressure(root.xylem.state.vc, config.KR_THRESHOLD)
        root.xylem.auxil.connected = false;

        # if at non-steady state, set the flow rate to be the sum of the buffer system so that flow from the soil is zero
        if root.xylem.auxil isa XylemHydraulicsAuxilNSS
            sol =  sum(root.xylem.auxil.flow_buffer);
        else
            sol = FT(0);
        end;
    else
        root.xylem.auxil.connected = true;

        # 1. set a max and min flow rate to use a bisection method to find the root flow rate
        p = abs(soil.auxil.ψ - junction.auxil.pressure - ρg_MPa(FT) * root.xylem.state.Δh);
        k = 1 / (1 / (root.rhizosphere.state.k_max * root.xylem.state.area) + 1 / (root.xylem.state.k_max * root.xylem.state.area / root.xylem.state.l));
        f_max = k * p;
        f_min = -f_max;

        # 2. set up a target function to find zero
        @inline root_pressure_diff(x::FT) where {FT} = (
            set_flow_profile!(root.xylem, x);
            root_pressure_profile!(soil, root, junction);

            return root.xylem.auxil.pressure[end] - junction.auxil.pressure
        );

        # 3. define method and solve for the root flow rate
        ms = NewtonBisectionMethod{FT}(x_min=f_min, x_max=f_max, x_ini=(f_min+f_max)/2);
        stol = SolutionTolerance{FT}(eps(FT)*100, 50);

        # 4. solve for the target flow rate
        sol = find_zero(root_pressure_diff, ms, stol);
    end;

    # 5. set up the flow rate using the solution from the solver
    set_flow_profile!(root.xylem, sol);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function root_flow_profiles!
#     2023-Sep-30: add root flow out into junction.auxil.∂w∂t
#
#######################################################################################################################################################################################################
"""

    root_flow_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}

Set up root flow profile for each root, given
- `config` `SPACConfiguration` type struct
- `spac` `BulkSPAC` type struct

"""
function root_flow_profiles!(config::SPACConfiguration{FT}, spac::BulkSPAC{FT}) where {FT}
    junction = spac.plant.junction;
    rindex = spac.plant.roots_index;
    roots = spac.plant.roots;
    soils = spac.soils;

    for i in eachindex(roots)
        root_flow_profile!(config, roots[i], soils[rindex[i]], junction);
        junction.auxil.∂w∂t += flow_out(roots[i]);
    end;

    return nothing
end;
