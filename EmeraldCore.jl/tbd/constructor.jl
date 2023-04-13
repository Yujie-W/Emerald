#=
TODO: move it to SoilPlantAirContinuum.jl
#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-May-25: add constructor function
#     2022-May-25: use Root and Stem structures with temperatures
#     2022-May-31: rename _qs to _fs
#     2022-May-31: add steady state mode option to input options
#     2022-Jun-15: fix documentation
#     2022-Jun-29: rename struct to MultiLayerSPAC, and use Leaves2D
#     2022-Jun-29: add CANOPY, Z, AIR, WLSET, LHA, ANGLES, SOIL, RAD_LW, RAD_SW, Φ_PHOTON to SPAC
#     2022-Jul-14: add area to constructor function
#
#######################################################################################################################################################################################################
"""

    MultiLayerSPAC{FT}(
                psm::String,
                area::Number = 100,
                wls::WaveLengthSet{FT} = WaveLengthSet{FT}();
                zs::Vector = [-1,6,12],
                zss::Vector = collect(0:-0.25:-2),
                zas::Vector = collect(0:0.5:13),
                ssm::Bool = true
    ) where {FT<:AbstractFloat}

Construct a SPAC system for monospecies tree system, given
- `psm` Photosynthesis model, must be C3, C4, or C3Cytochrome (note: there are C4 shrubs)
- `wls` [`WaveLengthSet`](@ref) type structure that determines the dimensions of leaf parameters
- `zs` Vector of Maximal root depth (negative value), trunk height, and canopy height
- `zss` Vector of soil layer boundaries starting from 0
- `zas` Vector of air layer boundaries starting from 0
- `ssm` Whether the flow rate is at steady state

---
# Examples
```julia
spac = MultiLayerSPAC{Float64}("C3");
```
"""
MultiLayerSPAC{FT}(
            psm::String,
            area::Number = 100,
            wls::WaveLengthSet{FT} = WaveLengthSet{FT}();
            zs::Vector = [-1,6,12],
            zss::Vector = collect(0:-0.25:-2),
            zas::Vector = collect(0:0.5:13),
            ssm::Bool = true
) where {FT<:AbstractFloat} = (
    @assert psm in ["C3", "C4", "C3Cytochrome"] "Photosynthesis model must be within [C3, C4, C3CytochromeModel]";

    # determine how many layers of roots
    _n_root = 0;
    _r_inds = Int[];
    for _i in eachindex(zss)
        if zss[_i] > zs[1]
            _n_root += 1;
            push!(_r_inds, _i);
        else
            break
        end;
    end;

    # determine how many layers of canopy
    _n_canopy= 0;
    _c_inds = Int[];
    for _i in 1:(length(zas)-1)
        # if the entire canopy is within the same layer
        if zas[_i] <= zs[2] < zs[3] <= zas[_i+1]
            _n_canopy += 1;
            push!(_c_inds, _i);
            break

        # if the lower canopy boundary (trunk top) is within the layer
        elseif zas[_i] < zs[2] < zas[_i+1] < zs[3]
            _n_canopy += 1;
            push!(_c_inds, _i);

        # if the layer is within the canopy
        elseif zs[2] <= zas[_i] < zas[_i+1] < zs[3]
            _n_canopy += 1;
            push!(_c_inds, _i);

        # if the upper canopy boundary is within the layer
        elseif zs[2] <= zas[_i] <= zs[3] < zas[_i+1]
            _n_canopy += 1;
            push!(_c_inds, _i-1);
            break

        # if the entire canopy is below the layer
        elseif zs[3] <= zas[_i] < zas[_i+1]
            break
        end;
    end;

    # create evenly distributed root system from top soil to deep soil
    _roots = Root{FT}[];
    for _i in _r_inds
        _Δh = abs(max(zss[_i+1], zs[1]) + zss[_i]) / 2;
        _hs = RootHydraulics{FT}(AREA = 1/_n_root, K_X = 25/_n_root, ΔH = _Δh);
        if !ssm _hs.FLOW = NonSteadyStateFlow{FT}(DIM_CAPACITY = 5) end;
        _rt = Root{FT}();
        _rt.HS = _hs;
        push!(_roots, _rt);
    end;

    # create trunk
    _hs = StemHydraulics{FT}(L = zs[2], ΔH = zs[2]);
    if !ssm
        _hs.FLOW = NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
    end;
    _trunk = Stem{FT}();
    _trunk.HS = _hs;

    # create branches from bottom canopy (trunk) to top canopy
    _branches = Stem{FT}[];
    for _i in _c_inds
        _Δh = (max(zas[_i], zs[2]) + min(zas[_i+1], zs[3])) / 2 - zs[2];
        _hs = StemHydraulics{FT}(AREA = 1/_n_canopy, ΔH = _Δh);
        if !ssm
            _hs.FLOW = NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
        end;
        _st = Stem{FT}();
        _st.HS = _hs;
        push!(_branches, _st);
    end;

    # create leaves from bottom canopy (trunk) to top canopy
    _leaves = [Leaves2D{FT}(psm, wls) for _i in 1:_n_canopy];
    for _leaf in _leaves
        _leaf.HS.AREA = 1500 / _n_canopy;
    end;

    # create canopy to use for radiative transfer
    _canopy = HyperspectralMLCanopy{FT}();

    # create air layers for all provided layers from bottom to top
    _airs = [AirLayer{FT}() for _i in 1:length(zas)-1];

    # create leaf hyperspectral absorption features
    _lha = HyperspectralAbsorption{FT}();

    # create sun sensor geometry
    _angles = SunSensorGeometry{FT}();

    # create soil
    _soil = Soil{FT}(ZS = zss);

    # create shortwave radiation
    _rad_sw = HyperspectralRadiation{FT}();

    # return plant
    return MultiLayerSPAC{FT}(
                _airs,              # AIR
                _angles,            # ANGLES
                _branches,          # BRANCHES
                _canopy,            # CANOPY
                _leaves,            # LEAVES
                _c_inds,            # LEAVES_INDEX
                _lha,               # LHA
                Meteorology{FT}(),  # METEO
                _n_canopy,          # N_CANOPY
                _n_root,            # N_ROOT
                100,                # RAD_LW
                _rad_sw,            # RAD_SW
                _roots,             # ROOTS
                _r_inds,            # ROOTS_INDEX
                _soil,              # SOIL
                _trunk,             # TRUNK
                wls,                # WLSET
                FT.(zs),            # Z
                FT.(zas),           # Z_AIR
                true,               # Φ_PHOTON
                zeros(FT,_n_root),  # _fs
                zeros(FT,_n_root),  # _ks
                zeros(FT,_n_root)   # _ps
    )
);
=#
