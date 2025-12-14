# This script is meant to speed up the SPAC model

using Emerald
using Revise

FT = Float64;

config = EmeraldLand.Namespace.SPACConfiguration(FT);




# default setting without binning PPAR values
spac = EmeraldLand.Namespace.BulkSPAC(config);
EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, 3600);
@info "GPP and SIF" EmeraldLand.SPAC.GPP(spac) EmeraldLand.SPAC.TROPOMI_SIF740(config, spac);

@time EmeraldLand.SPAC.initialize_spac!(config, spac);
@time EmeraldLand.SPAC.spac!(config, spac, 3600);




# the case using qL based fluorescence model
spac = EmeraldLand.Namespace.BulkSPAC(config);
for l in spac.plant.leaves
    l.photosystem.trait.FLM = EmeraldLand.Namespace.QLFluorescenceModelC3(FT);
end;
EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, 3600);
@info "GPP and SIF" EmeraldLand.SPAC.GPP(spac) EmeraldLand.SPAC.TROPOMI_SIF740(config, spac);

@time EmeraldLand.SPAC.initialize_spac!(config, spac);
@time EmeraldLand.SPAC.spac!(config, spac, 3600);




# the case using qL based fluorescence model Han
spac = EmeraldLand.Namespace.BulkSPAC(config);
for l in spac.plant.leaves
    l.photosystem.trait.FLM = EmeraldLand.Namespace.QLFluorescenceModelHanC3(FT);
end;
EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, 3600);
@info "GPP and SIF" EmeraldLand.SPAC.GPP(spac) EmeraldLand.SPAC.TROPOMI_SIF740(config, spac);

@time EmeraldLand.SPAC.initialize_spac!(config, spac);
@time EmeraldLand.SPAC.spac!(config, spac, 3600);




# the case using Cytochrome C3 model
spac = EmeraldLand.Namespace.BulkSPAC(config);
for l in spac.plant.leaves
    l.photosystem.trait.AJM = EmeraldLand.Namespace.AjMethodC3VqmaxPi();
    l.photosystem.trait.COLIMIT_J = EmeraldLand.Namespace.SerialColimit{FT}();
    l.photosystem.trait.FLM = EmeraldLand.Namespace.CytochromeFluorescenceModel{FT}();
end;
EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, 3600);
@info "GPP and SIF" EmeraldLand.SPAC.GPP(spac) EmeraldLand.SPAC.TROPOMI_SIF740(config, spac);

@time EmeraldLand.SPAC.initialize_spac!(config, spac);
@time EmeraldLand.SPAC.spac!(config, spac, 3600);




# the case using C4CLM model
spac = EmeraldLand.Namespace.BulkSPAC(config);
for l in spac.plant.leaves
    l.photosystem.trait = EmeraldLand.Namespace.GeneralC4Trait{FT}();
    l.photosystem.state = EmeraldLand.Namespace.C4State{FT}();
end;
EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, 3600);
@info "GPP and SIF" EmeraldLand.SPAC.GPP(spac) EmeraldLand.SPAC.TROPOMI_SIF740(config, spac);

@time EmeraldLand.SPAC.initialize_spac!(config, spac);
@time EmeraldLand.SPAC.spac!(config, spac, 3600);




# the case using C4VJP model
spac = EmeraldLand.Namespace.BulkSPAC(config);
for l in spac.plant.leaves
    l.photosystem.trait = EmeraldLand.Namespace.GeneralC4Trait{FT}();
    l.photosystem.trait.APM = EmeraldLand.Namespace.ApMethodC4VpmaxPi();
    l.photosystem.state = EmeraldLand.Namespace.C4State{FT}();
end;
EmeraldLand.SPAC.initialize_spac!(config, spac);
EmeraldLand.SPAC.spac!(config, spac, 3600);
@info "GPP and SIF" EmeraldLand.SPAC.GPP(spac) EmeraldLand.SPAC.TROPOMI_SIF740(config, spac);

@time EmeraldLand.SPAC.initialize_spac!(config, spac);
@time EmeraldLand.SPAC.spac!(config, spac, 3600);




# setting with binning PPAR values
config_b = EmeraldLand.Namespace.SPACConfiguration(FT);
config_b.DIM_PPAR_BINS = 20;
spac_b = EmeraldLand.Namespace.BulkSPAC(config_b);
EmeraldLand.SPAC.initialize_spac!(config_b, spac_b);
EmeraldLand.SPAC.spac!(config_b, spac_b, 3600);
@info "GPP and SIF" EmeraldLand.SPAC.GPP(spac_b) EmeraldLand.SPAC.TROPOMI_SIF740(config_b, spac_b);

@time EmeraldLand.SPAC.initialize_spac!(config_b, spac_b);
@time EmeraldLand.SPAC.spac!(config_b, spac_b, 3600);
