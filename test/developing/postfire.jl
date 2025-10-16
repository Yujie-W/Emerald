using Emerald
using Revise


FT = Float64;

config_1 = EmeraldLand.Namespace.SPACConfiguration(FT);
config_2 = EmeraldLand.Namespace.SPACConfiguration(FT);
spac_1 = EmeraldLand.Namespace.BulkSPAC(config);
spac_2 = EmeraldLand.Namespace.BulkSPAC(config);

for i in eachindex(spac_1.airs)
    spac_1.airs[i].auxil.wind = 1;
    spac_2.airs[i].auxil.wind = 1;
end;

EmeraldLand.SPAC.prescribe_traits!(config_1, spac_1; lai=0.5, sai=0.01);
EmeraldLand.SPAC.prescribe_traits!(config_1, spac_2; lai=0.5, sai=0.01);

config_2.SOIL_ALBEDO = Emerald.EmeraldLand.Namespace.SoilAlbedoPrescribe();

EmeraldLand.SPAC.initialize_spac!(config_1, spac_1);
EmeraldLand.SPAC.initialize_spac!(config_1, spac_2);

EmeraldLand.SPAC.spac!(config_1, spac_1, FT(3600));
EmeraldLand.SPAC.spac!(config_1, spac_2, FT(3600));
@info "comparison" spac_1.canopy.structure.auxil.lwꜛ[1] spac_2.canopy.structure.auxil.lwꜛ[1];

spac_2.soil_bulk.auxil.ρ_sw[config_2.SPECTRA.IΛ_PAR] .*= 0.1;



begin
    for i in 1:100
        EmeraldLand.SPAC.spac!(config_2, spac_1, FT(600));
        EmeraldLand.SPAC.spac!(config_2, spac_2, FT(600));
    end;
    println("LW out per layer");
    nlayer = length(spac_1.plant.leaves);
    for i in 1:nlayer
        @printf("%.1f    %.1f    %.1f    %.1f\n", spac_1.canopy.structure.auxil.lwꜛ[i], spac_2.canopy.structure.auxil.lwꜛ[i], spac_1.plant.leaves[nlayer+1-i].energy.s_aux.t, spac_2.plant.leaves[nlayer+1-i].energy.s_aux.t);
    end;
    @printf("%.1f    %.1f    %.1f    %.1f\n", spac_1.canopy.structure.auxil.lwꜛ[end], spac_2.canopy.structure.auxil.lwꜛ[end], spac_1.soils[1].s_aux.t, spac_2.soils[1].s_aux.t);
end;















begin
    config_debug = EmeraldLand.Namespace.SPACConfiguration(FT);
    spac_debug = EmeraldLand.Namespace.BulkSPAC(config);

    for i in eachindex(spac_1.airs)
        spac_debug.airs[i].auxil.wind = 0.1;
    end;

    EmeraldLand.SPAC.prescribe_traits!(config_debug, spac_debug; lai=2.0, sai=0.05);
    EmeraldLand.SPAC.initialize_spac!(config_debug, spac_debug);
    EmeraldLand.SPAC.spac!(config_debug, spac_debug, FT(1));

    EmeraldLand.SPAC.spac!(config_debug, spac_debug, FT(1));
    println("LW out per layer");
    nlayer = length(spac_debug.plant.leaves);
    for i in 1:nlayer
        @printf("%.1f    %.1f\n", spac_debug.canopy.structure.auxil.lwꜛ[i], spac_debug.plant.leaves[nlayer+1-i].energy.s_aux.t);
    end;
end;
