import Emerald.SPAC as ESPAC
import Emerald.Namespace as ENS


FT = Float64;
config = ENS.SPACConfig(FT; dataset = ENS.OLD_PHI_2021_1NM);
spac = ENS.BulkSPAC(config; air_bounds=collect(0:0.1:13));
ESPAC.prescribe_traits!(config, spac; lai=4.0, ci=1, sai=0);
ESPAC.initialize_spac!(config, spac);
ESPAC.spac!(config, spac, FT(360));


# compare the SIF at the first layer
sen_geo = spac.canopy.sensor_geometry;
sun_geo = spac.canopy.sun_geometry;
begin
    @show (sun_geo.auxil.e_sifꜜ_layer[:,1] .+ sun_geo.auxil.e_sifꜛ_layer[:,1]) ./ sun_geo.auxil.e_sif_chl[:,1];
end;
