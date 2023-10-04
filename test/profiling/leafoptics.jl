import Emerald.EmeraldLand.LeafOptics as LO
import Emerald.EmeraldLand.Namespace as NS

using Profile

config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
bio = NS.LeafBio(config);
LO.leaf_spectra!(config, bio, 5.0, 40.0; N = 10);
Profile.clear()

LO.leaf_spectra!(config, bio, 5.0, 40.0; N = 10);

exit()
