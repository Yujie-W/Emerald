using Emerald;
using Test;


@testset "Change SPAC PAR Definition" begin
    FT = Float64;

    # By default, the PAR definition is [300,750] after we include the UV and NIR bands into account.
    config_uvparfr = EmeraldLand.Namespace.SPACConfiguration{FT}(SPECTRA = EmeraldLand.Namespace.ReferenceSpectra{FT}());
    spac_uvparfr = EmeraldLand.Namespace.BulkSPAC(config_uvparfr);
    EmeraldLand.SPAC.initialize_spac!(config_uvparfr, spac_uvparfr);
    EmeraldLand.SPAC.spac!(config_uvparfr, spac_uvparfr, FT(600));
    @test true;

    # However, users can still use the old PAR definition [400,700] by setting the following
    # Note here that WL_PAR is the definition of PAR range, and WL_PAR_700 need to be changed as well to avoid any problem.
    config_paronly = EmeraldLand.Namespace.SPACConfiguration{FT}(SPECTRA = EmeraldLand.Namespace.ReferenceSpectra{FT}(WL_PAR = FT[400,700], WL_PAR_700 = [400,700]));
    spac_paronly = EmeraldLand.Namespace.BulkSPAC(config_paronly);
    EmeraldLand.SPAC.initialize_spac!(config_paronly, spac_paronly);
    EmeraldLand.SPAC.spac!(config_paronly, spac_paronly, FT(600));
    @test true;
end;
