#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2023-Mar-25: add a simple struct as in EmeraldEarth.LandDatasets
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure to save labels from GriddingMachine (year must be provided)

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct GriddingMachineLabels
    "Which year do the datasets apply (when applicable)"
    year::Int
    "GriddingMachine.jl tag for soil color class"
    tag_s_cc::String = "SC_2X_1Y_V1"
    "GriddingMachine.jl tag for soil van Genuchten parameters"
    tag_s_α::String = "SOIL_VGA_12X_1Y_V1"
    "GriddingMachine.jl tag for soil van Genuchten parameters"
    tag_s_n::String = "SOIL_VGN_12X_1Y_V1"
    "GriddingMachine.jl tag for soil van Genuchten parameters"
    tag_s_Θr::String = "SOIL_SWCR_12X_1Y_V1"
    "GriddingMachine.jl tag for soil van Genuchten parameters"
    tag_s_Θs::String = "SOIL_SWCS_12X_1Y_V1"
    "GriddingMachine.jl tag for canopy height"
    tag_p_ch::String = "CH_20X_1Y_V1"
    "GriddingMachine.jl tag for chlorophyll content"
    tag_p_chl::String = "CHL_2X_7D_V1"
    "GriddingMachine.jl tag for clumping index"
    tag_p_ci::String = "CI_2X_1M_V3"
    "GriddingMachine.jl tag for leaf area index"
    tag_p_lai::String = "LAI_MODIS_2X_8D_$(year)_V1"
    "GriddingMachine.jl tag for specific leaf area"
    tag_p_sla::String = "SLA_2X_1Y_V1"
    "GriddingMachine.jl tag for Vcmax"
    tag_p_vcm::String = "VCMAX_2X_1Y_V2"
    "GriddingMachine.jl tag for elevation"
    tag_t_ele::String = "ELEV_4X_1Y_V1"
    "GriddingMachine.jl tag for land mask"
    tag_t_lm::String = "LM_4X_1Y_V1"
    "GriddingMachine.jl tag for PFT"
    tag_t_pft::String = "PFT_2X_1Y_V1"
end


# CLM5 settings
CLM5_PFTG = [0, 2.35, 2.35, 2.35, 4.12, 4.12, 4.45, 4.45, 4.45, 4.7, 4.7, 4.7, 2.22, 5.25, 1.62, 5.79, 5.79] .* sqrt(1000);
CLM5_PFTS = ["not_vegetated",
             "needleleaf_evergreen_temperate",
             "needleleaf_evergreen_boreal",
             "needleleaf_deciduous_boreal",
             "broadleaf_evergreen_tropical",
             "broadleaf_evergreen_temperate",
             "broadleaf_deciduous_tropical",
             "broadleaf_deciduous_temperate",
             "broadleaf_deciduous_boreal",
             "evergreen_shrub",
             "deciduous_temperate_shrub",
             "deciduous_boreal_shrub",
             "c3_arctic_grass",
             "c3_non-arctic_grass",
             "c4_grass",
             "c3_crop",
             "c3_irrigated"];


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-25: add function to get dict so as to create a spac
#     2023-Mar-25: add three more fields to use ERA weather driver
#
#######################################################################################################################################################################################################
"""

    gm_dict(dts::GriddingMachineLabels, lat::Number, lon::Number)

Create a dict that contains data from GriddingMachine, given
- `dts` `GriddingMachineLabels` for GriddingMachine labels
- `lat` Latitude
- `lon` Longitude

"""
function gm_dict(dts::GriddingMachineLabels, lat::Number, lon::Number)
    # reading data from GriddingMachine
    if (read_LUT(query_collection(dts.tag_t_lm), lat, lon)[1] > 0)
        _lais = read_LUT(query_collection(dts.tag_p_lai), lat, lon)[1];
        if (nanmax(_lais) > 0)
            _ind_c3 = [2:14;16;17];
            _ccs = read_csv("$(@__DIR__)/../data/CO2-1Y.csv");
            _co2 = _ccs.MEAN[findfirst(_ccs.YEAR .== dts.year)];
            _pfts = read_LUT(query_collection(dts.tag_t_pft), lat, lon)[1];
            _g = CLM5_PFTG[_ind_c3]' * _pfts[_ind_c3] / sum(_pfts[_ind_c3]);
            _g1 = isnan(_g) ? nanmean(CLM5_PFTG[_ind_c3]) : _g;

            return Dict{String,Any}(
                        "C3C4"          => "C3",
                        "CANOPY_HEIGHT" => read_LUT(query_collection(dts.tag_p_ch), lat, lon)[1],
                        "CHLOROPHYLL"   => read_LUT(query_collection(dts.tag_p_chl), lat, lon)[1],
                        "CLUMPING"      => read_LUT(query_collection(dts.tag_p_ci), lat, lon)[1],
                        "CO2"           => _co2,
                        "ELEVATION"     => read_LUT(query_collection(dts.tag_t_ele), lat, lon)[1],
                        "FT"            => Float64,
                        "LAI"           => _lais,
                        "LAT_INDEX"     => lat_ind(lat; res = 1),
                        "LATITUDE"      => lat,
                        "LMA"           => 1 / 10 / read_LUT(query_collection(dts.tag_p_sla), lat, lon)[1],
                        "LON_INDEX"     => lon_ind(lon; res = 1),
                        "LONGITUDE"     => lon,
                        "MEDLYN_G1"     => _g1,
                        "RESO_SPACE"    => "1X",
                        "SOIL_COLOR"    => min(20, max(1, Int(floor(read_LUT(query_collection(dts.tag_s_cc), lat, lon)[1])))),
                        "SOIL_N"        => read_LUT(query_collection(dts.tag_s_n), lat, lon)[1],
                        "SOIL_α"        => read_LUT(query_collection(dts.tag_s_α), lat, lon)[1],
                        "SOIL_ΘR"       => read_LUT(query_collection(dts.tag_s_Θr), lat, lon)[1],
                        "SOIL_ΘS"       => read_LUT(query_collection(dts.tag_s_Θs), lat, lon)[1],
                        "VCMAX25"       => read_LUT(query_collection(dts.tag_p_vcm), lat, lon)[1],
                        "YEAR"          => dts.year,
            );
        end;
    end;

    return error("Not Land or LAI = 0...")
end
