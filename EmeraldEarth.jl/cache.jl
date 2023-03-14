#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the CACHE_SPAC
#     2023-Mar-13: initialize CACHE_STATE at the same time
#
#######################################################################################################################################################################################################
"""

    initialize_cache!(FT)

Initialize the global parameter `CACHE_SPAC`, given
- `FT` Floating type

"""
function initialize_cache!(FT)
    global CACHE_SPAC, CACHE_STATE;

    # create a SPAC to work on
    _z_canopy = FT(10);
    CACHE_SPAC = MonoMLTreeSPAC{FT}(
                DIM_AIR      = 25,
                DIM_LAYER    = 10,
                DIM_ROOT     = 4,
                LATITUDE     = 0,
                LONGITUDE    = 0,
                LEAVES_INDEX = collect(11:20),
                ROOTS_INDEX  = collect(1:4),
                Z            = [-2, _z_canopy/2, _z_canopy],
                Z_AIR        = collect(0:21) * _z_canopy / 20,
                SOIL         = Soil{FT}(DIM_SOIL = 4, ZS = [0, -0.1, -0.35, -1, -3]));

    # set hydraulic traits to very high so as to not triggering NaN (they do not impact result anyway)
    for _organ in [CACHE_SPAC.LEAVES; CACHE_SPAC.BRANCHES; CACHE_SPAC.TRUNK; CACHE_SPAC.ROOTS]
        _organ.HS.VC.B = 3;
        _organ.HS.VC.C = 1;
    end;

    # update leaf mass per area and stomtal model
    @inline linear_p_soil(x) = min(1, max(eps(FT), 1 + x / 5));
    _bt = BetaFunction{FT}(FUNC = linear_p_soil, PARAM_X = BetaParameterPsoil(), PARAM_Y = BetaParameterG1());
    for _leaves in CACHE_SPAC.LEAVES
        _leaves.SM = MedlynSM{FT}(G0 = 0.005, β = _bt);
    end;

    # initialize the spac
    initialize!(CACHE_SPAC);

    # create a state struct based on the spac
    CACHE_STATE = MonoMLTreeSPACState{FT}(CACHE_SPAC);

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-13: add function to initialize the CACHE_SPAC
#     2023-Mar-13: add step to synchronize state variables into CACHE_SPAC
#
#######################################################################################################################################################################################################
"""

    synchronize_cache!(gm_params::Dict{String,Any}, wd_params::Dict{String,Any}, state::Union{Nothing,MonoMLTreeSPACState{FT}}) where {FT<:AbstractFloat}

Synchronize SPAC parameters from,
- `gm_params` Dict for GriddingMachine parameters
- `wd_params` Dict for weather drivers
- `state` `MonoMLTreeSPACState` for all state variables, or nothing

"""
function synchronize_cache!(gm_params::Dict{String,Any}, wd_params::Dict{String,Any}, state::Union{Nothing,MonoMLTreeSPACState})
    FT = gm_params["FT"];
    _z_canopy = max(FT(0.1), gm_params["CANOPY_HEIGHT"]);

    #
    # TODO: update canopy height for plant hydraulic system based on _z_canopy
    #

    # update the values in the CACHE_SPAC, use .= for arrays
    global CACHE_SPAC;
    CACHE_SPAC.LATITUDE = gm_params["LATITUDE"];
    CACHE_SPAC.LONGITUDE = gm_params["LONGITUDE"];
    CACHE_SPAC.ELEVATION = gm_params["ELEVATION"];
    CACHE_SPAC.Z .= [-2, _z_canopy/2, _z_canopy];
    CACHE_SPAC.Z_AIR .= collect(0:21) * _z_canopy / 20;
    CACHE_SPAC.SOIL.COLOR = gm_params["SOIL_COLOR"];

    # update soil type information per layer
    for _i in eachindex(CACHE_SPAC.SOIL.LAYERS)
        # TODO: add a line to parameterize K_MAX
        CACHE_SPAC.SOIL.LAYERS[_i].VC.α = gm_params["SOIL_α"][_i];
        CACHE_SPAC.SOIL.LAYERS[_i].VC.N = gm_params["SOIL_N"][_i];
        CACHE_SPAC.SOIL.LAYERS[_i].VC.M = 1 - 1 / CACHE_SPAC.SOIL.LAYERS[_i].VC.N;
        CACHE_SPAC.SOIL.LAYERS[_i].VC.Θ_RES = gm_params["SOIL_ΘR"][_i];
        CACHE_SPAC.SOIL.LAYERS[_i].VC.Θ_SAT = gm_params["SOIL_ΘS"][_i];
    end;

    # update leaf mass per area and stomtal model
    for _leaves in CACHE_SPAC.LEAVES
        _leaves.BIO.lma = gm_params["LMA"];
        _leaves.SM.G1 = gm_params["MEDLYN_G1"];
    end;

    # update environmental conditions
    for _alayer in CACHE_SPAC.AIR
        _alayer.P_AIR = wd_params["P_ATM"];
        update!(_alayer; t = wd_params["T_AIR"], vpd = wd_params["VPD"], wind = wd_params["WIND"]);
    end;

    # sync the environmental conditions per layer for CO₂ concentration
    if !isnothing(gm_params["CO2"])
        for _alayer in CACHE_SPAC.AIR
            update!(_alayer; f_CO₂ = gm_params["CO2"]);
        end;
    end;

    # update shortwave radiation
    _in_dir = CACHE_SPAC.RAD_SW_REF.e_direct' * CACHE_SPAC.CANOPY.WLSET.ΔΛ / 1000;
    _in_dif = CACHE_SPAC.RAD_SW_REF.e_diffuse' * CACHE_SPAC.CANOPY.WLSET.ΔΛ / 1000;
    CACHE_SPAC.RAD_SW.e_direct  .= CACHE_SPAC.RAD_SW_REF.e_direct  .* wd_params["RAD_DIR"] ./ _in_dir;
    CACHE_SPAC.RAD_SW.e_diffuse .= CACHE_SPAC.RAD_SW_REF.e_diffuse .* wd_params["RAD_DIF"] ./ _in_dif;

    # update solar zenith angle based on the time
    CACHE_SPAC.ANGLES.sza = solar_zenith_angle(CACHE_SPAC.LATITUDE, FT(wd_params["FDOY"]));

    # prescribe soil water content
    if "SWC" in keys(wd_params)
        update!(CACHE_SPAC; swcs = wd_params["SWC"], t_soils = wd_params["T_SOIL"]);
    end;

    # synchronize the state if state is not nothing, otherwise set all values to NaN (do thing before prescribing T_SKIN)
    if !isnothing(state)
        spac_state!(state, CACHE_SPAC);
    else
        CACHE_SPAC.MEMORY.tem .= NaN;
    end;

    # prescribe leaf temperature from skin temperature
    # TODO: add CACHE_SPAC.MEMORY.tem as prognostic variable
    if "T_SKIN" in keys(wd_params)
        push!(CACHE_SPAC.MEMORY.tem, wd_params["T_SKIN"]);
        if length(CACHE_SPAC.MEMORY.tem) > 240 deleteat!(CACHE_SPAC.MEMORY.tem,1) end;
        update!(CACHE_SPAC; t_leaf = wd_params["T_SKIN"], t_clm = nanmean(CACHE_SPAC.MEMORY.tem));
    end;

    # synchronize LAI, CHL, and CI
    _iday = Int(floor(wd_params["INDEX"] / 24)) + 1;
    _chl = griddingmachine_data(gm_params["CHLOROPHYLL"], gm_params["YEAR"], _iday);
    _cli = griddingmachine_data(gm_params["CLUMPING"], gm_params["YEAR"], _iday);
    _lai = griddingmachine_data(gm_params["LAI"], gm_params["YEAR"], _iday);
    _vcm = griddingmachine_data(gm_params["VCMAX25"], gm_params["YEAR"], _iday);

    # update clumping index, LAI, Vcmax, and Chl
    CACHE_SPAC.CANOPY.ci = _cli;
    CACHE_SPAC.CANOPY.Ω_A = _cli;
    update!(CACHE_SPAC; cab = _chl, car = _chl / 7, lai =_lai, vcmax = _vcm, vcmax_expo = 0.3);

    return nothing
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to determine the day bounds of input GriddingMachine drivers
#
#######################################################################################################################################################################################################
"""

    griddingmachine_data_index(n::Int, year::Int, d::Int)

Return the index of data, given
- `n` Number of GriddingMachine data time index
- `year` Year
- `d` Day number

"""
function griddingmachine_data(data::Vector, year::Int, d::Int)
    _bounds = [0,367];
    _n = length(data);
    if _n == 1
        _bounds = [0,367]
    elseif _n == 12
        _bounds = isleapyear(year) ? MDAYS_LEAP : MDAYS;
    elseif _n == 46
        _bounds = [collect(0:8:361); 367]
    elseif _n == 52
        _bounds = [collect(0:7:361); 367]
    else
        @error "This temporal resolution is not supported: $(_n)!";
    end

    _ind = findfirst(d .<= _bounds) - 1

    return data[_ind]
end
