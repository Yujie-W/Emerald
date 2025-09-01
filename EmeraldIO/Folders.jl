"""

This module is meant to set up the folder structure for Emerald.

"""
module Folders


# AmeriFlux
AMERIFLUX_FOLDER      = joinpath(homedir(), "DATASERVER/field/AmeriFlux");
AMERIFLUX_DATA        = joinpath(AMERIFLUX_FOLDER, "original");
AMERIFLUX_REPROCESSED = joinpath(AMERIFLUX_FOLDER, "reprocessed");

# ERA5
ERA5_SL_HOURLY = "/mnt/net/ormosia/DATASERVER/model/ERA5/SingleLevels/Hourly";

# FluxNet2015
FLUXNET2015_FOLDER      = joinpath(homedir(), "DATASERVER/field/Fluxnet2015");
FLUXNET2015_DATA        = joinpath(FLUXNET2015_FOLDER, "data");
FLUXNET2015_REPROCESSED = joinpath(FLUXNET2015_FOLDER, "reprocessed");

# Land
LAND_FOLDER = "/mnt/net/ormosia/DATASERVER/model/Emerald";
LAND_CACHE  = joinpath(LAND_FOLDER, "cache");
LAND_DRIVER = joinpath(LAND_FOLDER, "drivers");
LAND_RESULT = joinpath(LAND_FOLDER, "simulations");
LAND_SETUP  = joinpath(LAND_FOLDER, "setups");


end # module
