# This file contains the spectra information to use with the SPAC system
# By default, the reference spectra is to run at hyperspectral mode
# However, if is also possible to run the model at broadband mode, and we provide a constructor for this purpose (SIF will be disabled this way)

#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2021-Aug-04: refactor the structure with constants, variables, and temporary cache
#     2021-Aug-04: add concentrations and characteristic curves altogether
#     2021-Aug-10: add CBC and PRO supoort
#     2022-Jul-20: add field DATASET to struct
#     2023-Sep-11: add field ΔΛ_SIF
#     2023-Sep-13: add fields Φ_PSI and Φ_PSII
#     2023-Sep-19: add fields MAT_SOIL and SOLAR_RAD
#     2023-Oct-17: add field ρ_STEM
#     2024-Jan-23: change default minimum PAR to 300 nm
#     2024-Jan-23: add field for PAR using the 300-700 nm definition
#     2024-Jan-23: redefine IΛ_PAR_700 to be the index within IΛ_PAR so that we can simply subset the PPAR even though the definition of PAR is changed
#     2024-Feb-27: remove field DATASET and create a new constructor function that specifies the dataset location
#     2024-Nov-11: remove fields Λ_UPPER and Λ_LOWER
#     2024-Nov-11: add option to read only the selected wavelengths
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Immutable struct that contains leaf biophysical traits used to run leaf reflection and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ReferenceSpectra{FT<:AbstractFloat}
    # Wavelegnth bins and upper/lower boundaries
    "Wavelength (bins) `[nm]`"
    Λ::Vector{FT}
    "Differential wavelength `[nm]`"
    ΔΛ::Vector{FT}

    # Constant features for the leaf
    "Specific absorption coefficients of anthocyanin `[-]`"
    K_ANT::Vector{FT}
    "Specific absorption coefficients of senescent material (brown pigments) `[-]`"
    K_BROWN::Vector{FT}
    "Specific absorption coefficients of chlorophyll a and b `[-]`"
    K_CAB::Vector{FT}
    "Specific absorption coefficients of violaxanthin carotenoid `[-]`"
    K_CAR_V::Vector{FT}
    "Specific absorption coefficients of zeaxanthin carotenoid `[-]`"
    K_CAR_Z::Vector{FT}
    "Specific absorption coefficients of carbon-based constituents `[-]`"
    K_CBC::Vector{FT}
    "Specific absorption coefficients of water `[-]`"
    K_H₂O::Vector{FT}
    "Specific absorption coefficients of dry matter `[-]`"
    K_LMA::Vector{FT}
    "Specific absorption coefficients of protein `[-]`"
    K_PRO::Vector{FT}
    "Refractive index `[-]`"
    NR::Vector{FT}
    "Fluorescence yield of PS I and II probability function `[nm⁻¹]`"
    Φ_PS::Vector{FT}
    "Fluorescence yield of PS I probability function `[nm⁻¹]`"
    Φ_PSI::Vector{FT}
    "Fluorescence yield of PS II probability function `[nm⁻¹]`"
    Φ_PSII::Vector{FT}

    # Stem reflectance
    "Stem reflectance `[-]`"
    ρ_STEM::Vector{FT} = zeros(FT, length(Λ)) .+ 0.2

    # Variable features for the soil
    "A matrix of characteristic curves"
    MAT_SOIL::Matrix{FT}

    # Variable features for the radiation
    "Downwelling shortwave radiation reference spectrum"
    SOLAR_RAD::Matrix{FT}

    # Constants
    "Wavelength limits for NIR `[nm]`"
    WL_NIR::Vector{FT} = FT[700, 2500]
    "Wavelength limits for PAR `[nm]`"
    WL_PAR::Vector{FT} = FT[300, 750]
    "Wavelength limits for PAR using the 300-700 nm setup `[nm]`"
    WL_PAR_700::Vector{FT} = FT[300, 700]
    "Wavelength limits for SIF emission `[nm]`"
    WL_SIF::Vector{FT} = FT[640, 850]
    "Wavelength limits for SIF excitation `[nm]`"
    WL_SIFE::Vector{FT} = FT[300, 750]

    # Indices
    "Indicies of Λ_NIR in Λ"
    IΛ_NIR::Vector{Int} = findall( WL_NIR[1] .<= Λ .<= WL_NIR[2] )
    "Indicies of Λ_PAR in Λ"
    IΛ_PAR::Vector{Int} = findall( WL_PAR[1] .<= Λ .<= WL_PAR[2] )
    "Indicies of Λ_PAR_700 in Λ_PAR using the 300-700 nm definition"
    IΛ_PAR_700::Vector{Int} = findall( WL_PAR_700[1] .<= Λ[IΛ_PAR] .<= WL_PAR_700[2] )
    "Indicies of Λ_SIF in Λ"
    IΛ_SIF::Vector{Int} = findall( WL_SIF[1] .<= Λ .<= WL_SIF[2] )
    "Indicies of Λ_SIFE in Λ"
    IΛ_SIFE::Vector{Int} = findall( WL_SIFE[1] .<= Λ .<= WL_SIFE[2] )

    # Constants based on the ones above
    "Differential wavelength for PAR `[nm]`"
    ΔΛ_PAR::Vector{FT} = ΔΛ[IΛ_PAR]
    "Differential wavelength for PAR using 300-700 nm definition `[nm]`"
    ΔΛ_PAR_700::Vector{FT} = ΔΛ_PAR[IΛ_PAR_700]
    "Differential wavelength for SIF `[nm]`"
    ΔΛ_SIF::Vector{FT} = ΔΛ[IΛ_SIF]
    "Differential wavelength for SIF excitation `[nm]`"
    ΔΛ_SIFE::Vector{FT} = ΔΛ[IΛ_SIFE]
    "Wavelength bins for PAR `[nm]`"
    Λ_PAR::Vector{FT} = Λ[IΛ_PAR]
    "Wavelength bins for PAR using the 300-700 nm definition `[nm]`"
    Λ_PAR_700::Vector{FT} = Λ_PAR[IΛ_PAR_700]
    "Wavelength bins for SIF `[nm]`"
    Λ_SIF::Vector{FT} = Λ[IΛ_SIF]
    "Wavelength bins for SIF excitation `[nm]`"
    Λ_SIFE::Vector{FT} = Λ[IΛ_SIFE]
end;

"""
#
    ReferenceSpectra{FT}(jld2_file::String, dataset::String; wl_par::Vector = [300,750], wl_par_700::Vector = [300,700], wl_selection::Union{Nothing,Vector} = nothing) where {FT}

Create and return a ReferenceSpectra object, given
- `jld2_file` the JLD2 file name
- `dataset` the dataset name in the JLD2 file
- `wl_par` the wavelength range for PAR
- `wl_par_700` the wavelength range for PAR 700
- `wl_selection` the wavelength selection

"""
ReferenceSpectra{FT}(jld2_file::String, dataset::String; wl_par::Vector = [300,750], wl_par_700::Vector = [300,700], wl_selection::Union{Nothing,Vector} = nothing) where {FT} = (
    df = read_jld2(jld2_file, dataset);

    # if wl_selection is provided, loop through the wavelengths and find the target values and use ΔΛ = 1 nm
    if !isnothing(wl_selection)
        if !occursin("1NM", dataset) && !occursin("1nm", dataset)
            @warn "It is recommended to use 1 nm resolution data when specifying the wavelength selection!";
        end;

        wls            = FT.(sort(wl_selection));
        Λ_interp       = wls;
        ΔΛ_interp      = ones(FT, length(wls));
        K_ANT_interp   = similar(wls, FT);
        K_BROWN_interp = similar(wls, FT);
        K_CAB_interp   = similar(wls, FT);
        K_CAR_V_interp = similar(wls, FT);
        K_CAR_Z_interp = similar(wls, FT);
        K_CBC_interp   = similar(wls, FT);
        K_H₂O_interp   = similar(wls, FT);
        K_LMA_interp   = similar(wls, FT);
        K_PRO_interp   = similar(wls, FT);
        NR_interp      = similar(wls, FT);
        Φ_PS_interp    = similar(wls, FT);
        Φ_PSI_interp   = similar(wls, FT);
        Φ_PSII_interp  = similar(wls, FT);
        GSV_1_interp   = similar(wls, FT);
        GSV_2_interp   = similar(wls, FT);
        GSV_3_interp   = similar(wls, FT);
        GSV_4_interp   = similar(wls, FT);
        E_DIR_interp   = similar(wls, FT);
        E_DIFF_interp  = similar(wls, FT);

        # interpolate the data
        for i in eachindex(wls)
            K_ANT_interp[i] = interpolate_data(df.WL, df.K_ANT, wls[i]);
            K_BROWN_interp[i] = interpolate_data(df.WL, df.K_BROWN, wls[i]);
            K_CAB_interp[i] = interpolate_data(df.WL, df.K_CAB, wls[i]);
            K_CAR_V_interp[i] = interpolate_data(df.WL, df.K_CAR_V, wls[i]);
            K_CAR_Z_interp[i] = interpolate_data(df.WL, df.K_CAR_Z, wls[i]);
            K_CBC_interp[i] = interpolate_data(df.WL, df.K_CBC, wls[i]);
            K_H₂O_interp[i] = interpolate_data(df.WL, df.K_H₂O, wls[i]);
            K_LMA_interp[i] = interpolate_data(df.WL, df.K_LMA, wls[i]);
            K_PRO_interp[i] = interpolate_data(df.WL, df.K_PRO, wls[i]);
            NR_interp[i] = interpolate_data(df.WL, df.NR, wls[i]);
            Φ_PS_interp[i] = interpolate_data(df.WL, df.K_PS, wls[i]);
            Φ_PSI_interp[i] = interpolate_data(df.WL, df.K_PS1, wls[i]);
            Φ_PSII_interp[i] = interpolate_data(df.WL, df.K_PS2, wls[i]);
            GSV_1_interp[i] = interpolate_data(df.WL, df.GSV_1, wls[i]);
            GSV_2_interp[i] = interpolate_data(df.WL, df.GSV_2, wls[i]);
            GSV_3_interp[i] = interpolate_data(df.WL, df.GSV_3, wls[i]);
            GSV_4_interp[i] = interpolate_data(df.WL, df.GSV_4, wls[i]);
            E_DIR_interp[i] = interpolate_data(df.WL, df.E_DIR, wls[i]);
            E_DIFF_interp[i] = interpolate_data(df.WL, df.E_DIFF, wls[i]);
        end;

        # return the ReferenceSpectra object
        return ReferenceSpectra{FT}(
                    Λ          = Λ_interp,
                    ΔΛ         = ΔΛ_interp,
                    K_ANT      = K_ANT_interp,
                    K_BROWN    = K_BROWN_interp,
                    K_CAB      = K_CAB_interp,
                    K_CAR_V    = K_CAR_V_interp,
                    K_CAR_Z    = K_CAR_Z_interp,
                    K_CBC      = K_CBC_interp,
                    K_H₂O      = K_H₂O_interp,
                    K_LMA      = K_LMA_interp,
                    K_PRO      = K_PRO_interp,
                    NR         = NR_interp,
                    Φ_PS       = Φ_PS_interp,
                    Φ_PSI      = Φ_PSI_interp,
                    Φ_PSII     = Φ_PSII_interp,
                    MAT_SOIL   = FT[GSV_1_interp GSV_2_interp GSV_3_interp GSV_4_interp],
                    SOLAR_RAD  = FT[E_DIR_interp E_DIFF_interp],
                    WL_PAR     = wl_par,
                    WL_PAR_700 = wl_par_700,
        )
    end;

    # if wl_selection is not provided, use the default
    return ReferenceSpectra{FT}(
            Λ          = df.WL,
            ΔΛ         = df.WL_UPPER - df.WL_LOWER,
            K_ANT      = df.K_ANT,
            K_BROWN    = df.K_BROWN,
            K_CAB      = df.K_CAB,
            K_CAR_V    = df.K_CAR_V,
            K_CAR_Z    = df.K_CAR_Z,
            K_CBC      = df.K_CBC,
            K_H₂O      = df.K_H₂O,
            K_LMA      = df.K_LMA,
            K_PRO      = df.K_PRO,
            NR         = df.NR,
            Φ_PS       = df.K_PS,
            Φ_PSI      = df.K_PS1,
            Φ_PSII     = df.K_PS2,
            MAT_SOIL   = FT[df.GSV_1 df.GSV_2 df.GSV_3 df.GSV_4],
            SOLAR_RAD  = FT[df.E_DIR df.E_DIFF],
            WL_PAR     = wl_par,
            WL_PAR_700 = wl_par_700,
    )
);

"""

    broadband_spectra(FT)

Return a ReferenceSpectra object with broadband spectra, given
- `FT`: the floating point type to use

"""
broadband_spectra(FT) = ReferenceSpectra{FT}(
            Λ         = [550, 1600],
            Λ_LOWER   = [400, 700],
            Λ_UPPER   = [700, 2500],
            ΔΛ        = [300, 1800],
            K_ANT     = [0, 0],                             # will not be used
            K_BROWN   = [0, 0],                             # will not be used
            K_CAB     = [0, 0],                             # will not be used
            K_CAR_V   = [0, 0],                             # will not be used
            K_CAR_Z   = [0, 0],                             # will not be used
            K_CBC     = [0, 0],                             # will not be used
            K_H₂O     = [0, 0],                             # will not be used
            K_LMA     = [0, 0],                             # will not be used
            K_PRO     = [0, 0],                             # will not be used
            NR        = [1.5, 1.5],                         # will not be used
            Φ_PS      = [0, 0],                             # will not be used
            Φ_PSI     = [0, 0],                             # will not be used
            Φ_PSII    = [0, 0],                             # will not be used
            ρ_STEM    = [0.16, 0.39],
            MAT_SOIL  = [0 0 0 0],                          # will not be used
            SOLAR_RAD = [191.544 119.923; 277.752 55.076],
            WL_NIR    = [700, 2500],
            WL_PAR    = [400, 700],
            WL_SIF    = [640, 850],                         # will not be used
            WL_SIFE   = [400, 750],                         # will not be used
            IΛ_NIR    = [2],
            IΛ_PAR    = [1],
            IΛ_SIF    = [1],                                # will not be used
            IΛ_SIFE   = [1],                                # will not be used
            ΔΛ_PAR    = [300],
            ΔΛ_SIF    = [210],                              # will not be used
            ΔΛ_SIFE   = [350],                              # will not be used
            Λ_PAR     = [550],
            Λ_SIF     = [745],                              # will not be used
            Λ_SIFE    = [550],                              # will not be used
);
