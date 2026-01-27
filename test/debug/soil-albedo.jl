import Emerald.Namespace as ENS
import Emerald.SPAC as ESPAC


config = ENS.SPACConfig(Float64);
config.METHODS.SOIL_ALBEDO = ENS.SoilAlbedoPrescribe();
spac = ENS.BulkSPAC(config);

ESPAC.prescribe_traits!(config, spac; lai = 4.0, ci = 1.0, sai = 0, vcmax = 140, vertical_expo = 0.3, cab = 15);
for l in spac.plant.leaves
    l.bio.trait.lma = 0.012;
end;

exclusions = @. config.CONSTANTS.SPECTRA.Λ <= 500 || config.CONSTANTS.SPECTRA.Λ >= 525;
spac.meteo.rad_sw.e_dir[exclusions] .= 0.0;
spac.meteo.rad_sw.e_dif[exclusions] .= 0.0

ESPAC.initialize_spac!(config, spac);
ESPAC.spac!(config, spac, 360);
@show spac.soil_bulk.auxil.ρ_sw;
for l in spac.plant.leaves
    @show maximum(l.photosystem.auxil.ϕ_f);
end;
