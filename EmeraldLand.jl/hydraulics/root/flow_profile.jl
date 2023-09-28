# This function contains functions to set up root flow profile

#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-23: add function root_flow_profile!
#
#######################################################################################################################################################################################################
"""

    root_flow_profile!(root::Root{FT}, soil::SoilLayer{FT}, p_target::FT) where {FT}

Use solver to determine the root flow rate, given
- `root` `Root` type struct
- `soil` `SoilLayer` type struct
- `p_target` Target pressure at the end of root xylem

"""
function root_flow_profile!(root::Root{FT}, soil::SoilLayer{FT}, p_target::FT) where {FT}
    # 1. set a max and min flow rate to use a bisection method to find the root flow rate
    p = abs(soil.ψ - p_target - ρg_MPa(FT) * root.xylem.state.Δh);
    k = 1 / (1 / (root.rhizosphere.state.k_max * root.xylem.state.area) + 1 / (root.xylem.state.k_max * root.xylem.state.area / root.xylem.state.l));
    f_max = k * p;
    f_min = -f_max;

    # 2. set up a target function to find zero
    @inline root_pressure_diff(x::FT) where {FT} = (
        set_flow_profile!(root.xylem, x);
        root_pressure_profile!(root, soil);

        return root.xylem.auxil.pressure[end] - p_target
    );

    # 3. define method and solve for the root flow rate
    ms = NewtonBisectionMethod{FT}(x_min=f_min, x_max=f_max, x_ini=(f_min+f_max)/2);
    stol = SolutionTolerance{FT}(eps(FT)*100, 50);

    # 4. solve for the target flow rate
    sol = find_zero(root_pressure_diff, ms, stol);

    # 5. set up the flow rate using the solution from the solver
    set_flow_profile!(root.xylem, sol);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Sep-28: add function root_flow_profiles!
#
#######################################################################################################################################################################################################
"""

    root_flow_profiles!(spac::MultiLayerSPAC{FT}) where {FT}

Set up root flow profile for each root, given
- `spac` `MultiLayerSPAC` type struct

"""
function root_flow_profiles!(spac::MultiLayerSPAC{FT}) where {FT}
    (; JUNCTION, ROOTS, ROOTS_INDEX, SOIL) = spac;

    for i in eachindex(ROOTS)
        root_flow_profile!(ROOTS[i].NS, SOIL.LAYERS[ROOTS_INDEX[i]], JUNCTION.auxil.pressure);
    end;

    return nothing
end;
