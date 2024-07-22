using DataFrames;
using Emerald;
using LazyArtifacts;
using Test;


@testset "Generate Canopy SIF Spectrum" begin
    FT = Float64;

    # The default spectra using the old phi settings
    config = EmeraldLand.Namespace.SPACConfiguration(FT; dataset = EmeraldLand.Namespace.OLD_PHI_2021_1NM);
    config.Φ_SIF_WL = false;
    spac = EmeraldLand.Namespace.BulkSPAC(config);
    EmeraldLand.SPAC.initialize_spac!(config, spac);
    EmeraldLand.SPAC.spac!(config, spac, 3600);
    @test true;

    # Use the spectra used by Yinon Bar-On
    config_new = EmeraldLand.Namespace.SPACConfiguration(FT; dataset = EmeraldLand.Namespace.NEW_PHI_2021_1NM);
    config_new.Φ_SIF_WL = false;
    spac_new = EmeraldLand.Namespace.BulkSPAC(config_new);
    EmeraldLand.SPAC.initialize_spac!(config_new, spac_new);
    EmeraldLand.SPAC.spac!(config_new, spac_new, 3600);
    @test true;

    # Use my WL-dependent spectra
    config_def = EmeraldLand.Namespace.SPACConfiguration(FT; dataset = EmeraldLand.Namespace.OLD_PHI_2021_1NM);
    config_def.Φ_SIF_WL = true;
    spac_def = EmeraldLand.Namespace.BulkSPAC(config_def);
    EmeraldLand.SPAC.initialize_spac!(config_def, spac_def);
    EmeraldLand.SPAC.spac!(config_def, spac_def, 3600);
    @test true;

    # save the data as a CSV file
    #=
    df = DataFrame(
                WL = config.SPECTRA.Λ_SIF,
                SIF_OLD = spac.canopy.sensor_geometry.auxil.sif_obs,
                SIF_NEW = spac_new.canopy.sensor_geometry.auxil.sif_obs,
                SIF_DEF = spac_def.canopy.sensor_geometry.auxil.sif_obs);
    EmeraldIO.Text.save_csv!("sif-spectra.csv", df);
    =#
end;
