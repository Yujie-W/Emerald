using GriddingMachine.Indexer: LandDatasetLabels, grid_dict
using Emerald.Land: simulation!
using OrderedCollections: OrderedDict


simulation_yx!(settings::Union{Dict,OrderedDict}, lat::Number, lon::Number, year::Int; c3c4::String = "C3", saving::Union{Nothing,String} = nothing) =
(
    gmd = grid_dict(LandDatasetLabels(settings["GM_VERSION"], year), lat, lon);
    gmd["LAI"] .= 0.5;
    gmd["SAI"] .= 0.0;
    simulation!(settings, gmd; c3c4 = c3c4, saving = saving);
);

thread_func(p) = (
    settings, lat, lon, year, filename = p;
    simulation_yx!(settings, lat, lon, year; c3c4 = "C3", saving = filename);
);
