import Emerald.SPAC as ESPAC
import Emerald.Namespace as ENS


FT = Float64;
config = ENS.SPACConfig(FT; dataset = ENS.OLD_PHI_2021_1NM);
spac = ENS.BulkSPAC(config; air_bounds=collect(0:0.5:13));

tar_wl = 600;
spac.meteo.rad_sw.e_dir .*= 100;
spac.meteo.rad_sw.e_dif .*= 100;
spac.meteo.rad_sw.e_dir[config.CONSTANTS.SPECTRA.Λ .!= tar_wl] .= 0;
spac.meteo.rad_sw.e_dif[config.CONSTANTS.SPECTRA.Λ .!= tar_wl] .= 0;

ESPAC.prescribe_traits!(config, spac; lai=5.0, ci=1, sai=0);
ESPAC.initialize_spac!(config, spac);
ESPAC.spac!(config, spac, FT(360));


# compare the SIF at the first layer
sen_geo = spac.canopy.sensor_geometry;
sun_geo = spac.canopy.sun_geometry;
begin
    for i in eachindex(spac.plant.leaves)
        lratios = (sun_geo.auxil.e_sifꜜ_layer[:,i] .+ sun_geo.auxil.e_sifꜛ_layer[:,i]) ./ sun_geo.auxil.e_sif_chl[:,i];
        println(lratios[[50,100,150]]);
    end;
end;
