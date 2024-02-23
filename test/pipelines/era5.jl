# This file is meant to automatically download and regrid the ERA5 single levels data
using Emerald.EmeraldData.WeatherDrivers: fetch_ERA5_data!, regrid_ERA5!

# download ERA5 data
years = 2000:2022;
for year in years
    fetch_ERA5_data!(year);
end;

# regrid ERA5 data
for year in years
    regrid_ERA5!(year, 1);
end;
