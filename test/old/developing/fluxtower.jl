#
# this script is meant to test the new EmeraldData.FluxTower module
#
using Emerald;
using Emerald.EmeraldData.FluxTower: FluxTowerDataset, US_NR1, process_data!

# create a FluxTowerData struct
ftd = FluxTowerDataset(US_NR1());
dfout = process_data!(ftd; displaying = true, saving = false);
