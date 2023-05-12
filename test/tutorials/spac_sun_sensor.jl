@testset "Modify Sun Sensor Geometry" begin
    FT = Float64;
    config = EmeraldLand.Namespace.SPACConfiguration{FT}();
    spac = EmeraldLand.Namespace.MultiLayerSPAC{FT}();
    EmeraldLand.SPAC.initialize!(spac, config);
    EmeraldLand.SPAC.spac!(spac, config, FT(1));

    # Sun-sensor geometry impacts the signal that can be detected remotely.
    # There are several parameters related, and they are
    #     solar zenith angle,
    #     solar azimuth angle,
    #     viewer zenith angle, and
    #     viewer azimuth angle.
    # Users may customize these angles freely, and then run function spac! to compare the difference.
    spac.ANGLES.sza = 10;
    spac.ANGLES.saa = 180;
    spac.ANGLES.vza = 20;
    spac.ANGLES.vaa = 170;
    @test true;

    # We also have an embedded function to compute solar zenith angle based on latitude (e.g., 30°N) and solar time (day 1, 13:15; solar noon at 12:00).
    # For more details of this function, please refer to the documentation.
    spac.ANGLES.sza = EmeraldLand.EarthGeometry.solar_zenith_angle(30.0, 1.0, 13.0, 15.0);
    @test true;
end
