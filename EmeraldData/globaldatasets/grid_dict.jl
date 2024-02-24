#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-25: add function to get dict so as to create a spac
#     2023-Mar-25: add three more fields to use ERA weather driver
#     2023-Aug-25: interpolate the data after reading from GriddingMachine
#     2024-Feb-22: separate the grid_dict function from the gm_grids function that prepare all the grids
#     2024-Feb-23: use 0 for the plant-related fields for non-vegetated land
#     2024-Feb-23: set SAI to be 1/10 of the maximum LAI
#
#######################################################################################################################################################################################################
"""

    grid_dict(dts::LandDatasets{FT}, ilat::Int, ilon::Int; ccs::DataFrame = CCS) where {FT}
    grid_dict(dtl::LandDatasetLabels, year::Int, nx::Int, lat::Number, lon::Number; ccs::DataFrame = CCS)

Prepare a dictionary of GriddingMachine data to feed SPAC, given
- `dts` `LandDatasets` type data struct
- `ilat` latitude index
- `ilon` longitude index
- `ccs` CO2 concentration data (provided by default)
- `dtl` `LandDatasetLabels` type data struct
- `year` year of the datasets
- `nx` grid resolution (1/nx °)
- `lat` latitude
- `lon` longitude

"""
function grid_dict end;

grid_dict(dts::LandDatasets{FT}, ilat::Int, ilon::Int; ccs::DataFrame = CCS) where {FT} = (
    reso   = 1 / dts.LABELS.nx;
    co2    = ccs.MEAN[findfirst(ccs.YEAR .== dts.LABELS.year)];
    lmsk   = dts.t_lm[ilon,ilat,1];
    scolor = min(20, max(1, Int(floor(dts.s_cc[ilon,ilat,1]))));
    s_α    = dts.s_α[ilon,ilat,:];
    s_n    = dts.s_n[ilon,ilat,:];
    s_Θr   = dts.s_Θr[ilon,ilat,:];
    s_Θs   = dts.s_Θs[ilon,ilat,:];

    # return the grid dictionary if the grid is masked as soil
    if dts.mask_soil[ilon,ilat]
        return Dict{String,Any}(
                    "CANOPY_HEIGHT" => 0,
                    "CHLOROPHYLL"   => [0],
                    "CLUMPING"      => [1],
                    "CO2"           => co2,
                    "ELEVATION"     => dts.t_ele[ilon,ilat],
                    "FT"            => FT,
                    "G1_MEDLYN_C3"  => 0,
                    "G1_MEDLYN_C4"  => 0,
                    "LAI"           => [0],
                    "LAND_MASK"     => lmsk,
                    "LATITUDE"      => (ilat - 0.5) * reso - 90,
                    "LAT_INDEX"     => ilat,
                    "LMA"           => 0,
                    "LONGITUDE"     => (ilon - 0.5) * reso - 180,
                    "LON_INDEX"     => ilon,
                    "MESSAGE_LEVEL" => 0,
                    "PFT_FRACTIONS" => [0],
                    "RESO_SPACE"    => dts.LABELS.nx,
                    "SAI"           => 0,
                    "SOIL_COLOR"    => scolor,
                    "SOIL_N"        => s_n,
                    "SOIL_α"        => s_α,
                    "SOIL_ΘR"       => s_Θr,
                    "SOIL_ΘS"       => s_Θs,
                    "VCMAX25"       => 0,
                    "YEAR"          => dts.LABELS.year,
                    "ρ_NIR_C3"      => 0,
                    "ρ_NIR_C4"      => 0,
                    "ρ_PAR_C3"      => 0,
                    "ρ_PAR_C4"      => 0,
                    "τ_NIR_C3"      => 0,
                    "τ_NIR_C4"      => 0,
                    "τ_PAR_C3"      => 0,
                    "τ_PAR_C4"      => 0);
    end;

    # else return the grid dictionary if the grid is masked as plant
    chls  = dts.p_chl[ilon,ilat,:];
    cis   = dts.p_ci[ilon,ilat,:];
    lais  = dts.p_lai[ilon,ilat,:];
    lma   = 1 / dts.p_sla[ilon,ilat,1] / 10;
    pfts  = dts.t_pft[ilon,ilat,:];
    vcmax = dts.p_vcm[ilon,ilat,:];
    zc    = dts.p_ch[ilon,ilat,1];

    # gap fill the data for seasonal trends
    interpolate_data!(chls);
    interpolate_data!(cis);
    interpolate_data!(lais);
    interpolate_data!(vcmax);

    # compute g1 for Medlyn model
    ind_c3 = [2:14;16;17];
    ind_c4 = [15];

    g1_c3_medlyn = CLM5_PFTG[ind_c3]' * pfts[ind_c3] / sum(pfts[ind_c3]);
    g1_c4_medlyn = CLM5_PFTG[ind_c4]' * pfts[ind_c4] / sum(pfts[ind_c4]);
    if isnan(g1_c3_medlyn) g1_c3_medlyn = nanmean(CLM5_PFTG[ind_c3]) end;
    if isnan(g1_c4_medlyn) g1_c4_medlyn = nanmean(CLM5_PFTG[ind_c4]) end;

    ρ_par_c3 = CLM5_ρPAR[ind_c3]' * pfts[ind_c3] / sum(pfts[ind_c3]);
    τ_par_c3 = CLM5_τPAR[ind_c3]' * pfts[ind_c3] / sum(pfts[ind_c3]);
    ρ_nir_c3 = CLM5_ρNIR[ind_c3]' * pfts[ind_c3] / sum(pfts[ind_c3]);
    τ_nir_c3 = CLM5_τNIR[ind_c3]' * pfts[ind_c3] / sum(pfts[ind_c3]);
    ρ_par_c4 = CLM5_ρPAR[ind_c4]' * pfts[ind_c4] / sum(pfts[ind_c4]);
    τ_par_c4 = CLM5_τPAR[ind_c4]' * pfts[ind_c4] / sum(pfts[ind_c4]);
    ρ_nir_c4 = CLM5_ρNIR[ind_c4]' * pfts[ind_c4] / sum(pfts[ind_c4]);
    τ_nir_c4 = CLM5_τNIR[ind_c4]' * pfts[ind_c4] / sum(pfts[ind_c4]);

    return Dict{String,Any}(
                "CANOPY_HEIGHT" => max(0.1, zc),
                "CHLOROPHYLL"   => chls,
                "CLUMPING"      => cis,
                "CO2"           => co2,
                "ELEVATION"     => dts.t_ele[ilon,ilat],
                "FT"            => FT,
                "G1_MEDLYN_C3"  => g1_c3_medlyn,
                "G1_MEDLYN_C4"  => g1_c4_medlyn,
                "LAI"           => lais,
                "LAND_MASK"     => lmsk,
                "LATITUDE"      => (ilat - 0.5) * reso - 90,
                "LAT_INDEX"     => ilat,
                "LMA"           => lma,
                "LONGITUDE"     => (ilon - 0.5) * reso - 180,
                "LON_INDEX"     => ilon,
                "MESSAGE_LEVEL" => 0,
                "PFT_FRACTIONS" => pfts,
                "RESO_SPACE"    => dts.LABELS.nx,
                "SAI"           => nanmax(lais) / 10,
                "SOIL_COLOR"    => scolor,
                "SOIL_N"        => s_n,
                "SOIL_α"        => s_α,
                "SOIL_ΘR"       => s_Θr,
                "SOIL_ΘS"       => s_Θs,
                "VCMAX25"       => vcmax,
                "YEAR"          => dts.LABELS.year,
                "ρ_NIR_C3"      => ρ_nir_c3,
                "ρ_NIR_C4"      => ρ_nir_c4,
                "ρ_PAR_C3"      => ρ_par_c3,
                "ρ_PAR_C4"      => ρ_par_c4,
                "τ_NIR_C3"      => τ_nir_c3,
                "τ_NIR_C4"      => τ_nir_c4,
                "τ_PAR_C3"      => τ_par_c3,
                "τ_PAR_C4"      => τ_par_c4);
);

grid_dict(dtl::LandDatasetLabels, lat::Number, lon::Number; FT::DataType = Float64, ccs::DataFrame = CCS) = (
    lmsk = read_LUT(query_collection(dtl.tag_t_ele), lat, lon)[1];
    if !(lmsk > 0)
        return error("The target grid does not contain land!");
    end;

    lais = read_LUT(query_collection(dtl.tag_p_lai), lat, lon)[1];
    if !(nanmax(lais) > 0)
        return error("The target grid is not vegetated!");
    end;

    co2 = ccs.MEAN[findfirst(ccs.YEAR .== dtl.year)];
    scolor = min(20, max(1, Int(floor(read_LUT(query_collection(dtl.tag_s_cc), lat, lon)[1]))));
    s_α = read_LUT(query_collection(dtl.tag_s_n), lat, lon)[1];
    s_n = read_LUT(query_collection(dtl.tag_s_n), lat, lon)[1];
    s_Θr = read_LUT(query_collection(dtl.tag_s_Θr), lat, lon)[1];
    s_Θs = read_LUT(query_collection(dtl.tag_s_Θs), lat, lon)[1];

    # else return the grid dictionary if the grid is masked as plant
    chls = read_LUT(query_collection(dtl.tag_p_chl), lat, lon)[1];
    cis = read_LUT(query_collection(dtl.tag_p_ci), lat, lon)[1];
    lma = 1 / read_LUT(query_collection(dtl.tag_p_sla), lat, lon)[1] / 10;
    pfts = read_LUT(query_collection(dtl.tag_t_pft), lat, lon)[1];
    zc = read_LUT(query_collection(dtl.tag_p_ch), lat, lon)[1];

    if dtl.gm_tag == "gm3"
        vcmax = read_LUT(query_collection("VCMAX_2X_1Y_V2"), lat, lon)[1] .* 0.6;
    else
        vcmax = read_LUT(query_collection(dtl.tag_p_vcm), lat, lon)[1];
    end;

    # gap fill the data for seasonal trends
    interpolate_data!(chls);
    interpolate_data!(cis);
    interpolate_data!(lais);
    interpolate_data!(vcmax);

    # compute g1 for Medlyn model
    ind_c3 = [2:14;16;17];
    ind_c4 = [15];

    g1_c3_medlyn = CLM5_PFTG[ind_c3]' * pfts[ind_c3] / sum(pfts[ind_c3]);
    g1_c4_medlyn = CLM5_PFTG[ind_c4]' * pfts[ind_c4] / sum(pfts[ind_c4]);
    if isnan(g1_c3_medlyn) g1_c3_medlyn = nanmean(CLM5_PFTG[ind_c3]) end;
    if isnan(g1_c4_medlyn) g1_c4_medlyn = nanmean(CLM5_PFTG[ind_c4]) end;

    ρ_par_c3 = CLM5_ρPAR[ind_c3]' * pfts[ind_c3] / sum(pfts[ind_c3]);
    τ_par_c3 = CLM5_τPAR[ind_c3]' * pfts[ind_c3] / sum(pfts[ind_c3]);
    ρ_nir_c3 = CLM5_ρNIR[ind_c3]' * pfts[ind_c3] / sum(pfts[ind_c3]);
    τ_nir_c3 = CLM5_τNIR[ind_c3]' * pfts[ind_c3] / sum(pfts[ind_c3]);
    ρ_par_c4 = CLM5_ρPAR[ind_c4]' * pfts[ind_c4] / sum(pfts[ind_c4]);
    τ_par_c4 = CLM5_τPAR[ind_c4]' * pfts[ind_c4] / sum(pfts[ind_c4]);
    ρ_nir_c4 = CLM5_ρNIR[ind_c4]' * pfts[ind_c4] / sum(pfts[ind_c4]);
    τ_nir_c4 = CLM5_τNIR[ind_c4]' * pfts[ind_c4] / sum(pfts[ind_c4]);

    return Dict{String,Any}(
                "CANOPY_HEIGHT" => max(0.1, zc),
                "CHLOROPHYLL"   => chls,
                "CLUMPING"      => cis,
                "CO2"           => co2,
                "ELEVATION"     => read_LUT(query_collection(dtl.tag_t_ele), lat, lon)[1],
                "FT"            => FT,
                "G1_MEDLYN_C3"  => g1_c3_medlyn,
                "G1_MEDLYN_C4"  => g1_c4_medlyn,
                "LAI"           => lais,
                "LAND_MASK"     => lmsk,
                "LATITUDE"      => lat,
                "LAT_INDEX"     => lat_ind(lat; res = 1/dtl.nx),
                "LMA"           => lma,
                "LONGITUDE"     => lon,
                "LON_INDEX"     => lon_ind(lon; res = 1/dtl.nx),
                "MESSAGE_LEVEL" => 0,
                "RESO_SPACE"    => dtl.nx,
                "PFT_FRACTIONS" => pfts,
                "SAI"           => nanmax(lais) / 10,
                "SOIL_COLOR"    => scolor,
                "SOIL_N"        => s_n,
                "SOIL_α"        => s_α,
                "SOIL_ΘR"       => s_Θr,
                "SOIL_ΘS"       => s_Θs,
                "VCMAX25"       => vcmax,
                "YEAR"          => dtl.year,
                "ρ_NIR_C3"      => ρ_nir_c3,
                "ρ_NIR_C4"      => ρ_nir_c4,
                "ρ_PAR_C3"      => ρ_par_c3,
                "ρ_PAR_C4"      => ρ_par_c4,
                "τ_NIR_C3"      => τ_nir_c3,
                "τ_NIR_C4"      => τ_nir_c4,
                "τ_PAR_C3"      => τ_par_c3,
                "τ_PAR_C4"      => τ_par_c4);
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: migrate from research repo to Emerald
#     2023-Jun-15: add non-vegetated land in global simulations
#     2024-Feb-23: rename the function to grid_dict_mat after separating the grid_dict function
#
#######################################################################################################################################################################################################
"""

    grid_dict_mat(dts::LandDatasets{FT}) where {FT}

Prepare a matrix of GriddingMachine data to feed SPAC, given
- `dts` `LandDatasets` type data struct

"""
function grid_dict_mat(dts::LandDatasets{FT}) where {FT}
    # create a matrix of GriddingMachine data
    # TODO: add a step to verify the input datasets
    @tinfo "Preparing a matrix of GriddingMachine data to work on...";
    mat_gm = Matrix{Union{Nothing,Dict{String,Any}}}(nothing, size(dts.t_lm));
    for ilon in axes(dts.t_lm,1), ilat in axes(dts.t_lm,2)
        if dts.mask_spac[ilon,ilat] || dts.mask_soil[ilon,ilat]
            mat_gm[ilon,ilat] = grid_dict(dts, ilat, ilon);
        end;
    end;

    return mat_gm
end;
