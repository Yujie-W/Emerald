using Emerald
using Emerald.Namespace
using Emerald.SPAC


FT = Float64;
config = Emerald.Namespace.SPACConfig(FT);
#=
        config.METHODS.C3_AC_METHOD = Namespace.AcMethodC3VcmaxPi();
        config.METHODS.C3_AJ_METHOD = Namespace.AjMethodC3VqmaxPi();
        config.METHODS.C3_AP_METHOD = Namespace.ApMethodC3Vcmax();
        config.METHODS.COLIMIT_J = Namespace.SerialColimit();
        config.METHODS.TD_VCMAX_C3.ΔHA = 63000;
        config.METHODS.TD_VCMAX_C3.ΔHD = 204000;
        config.METHODS.TD_KQ = Namespace.ArrheniusPeak{Float64}(298.15, 300, 28500, 223500, 700);
        config.METHODS.TD_Γ.VAL_REF = 4.56;
        config.METHODS.TD_Γ.ΔHA = 11800;
        config.METHODS.TD_ηC.VAL_REF = 4 * 3 / 14;
        config.METHODS.TD_ηC.ΔHA = 28500;
        config.METHODS.TD_ηC.ΔHD = 223500;
        config.METHODS.TD_ηC.ΔSV = 700;
        config.METHODS.TD_ηL.VAL_REF = 3 * 3 / 14;
        config.METHODS.TD_ηL.ΔHA = 28500;
        config.METHODS.TD_ηL.ΔHD = 223500;
        config.METHODS.TD_ηL.ΔSV = 700;
config.METHODS.FLUORESCENCE_METHOD_C3 = Namespace.CytochromeFluorescenceModel();
=#

spac = Namespace.BulkSPAC(config);
SPAC.initialize_spac!(config, spac);
for i in 1:7200
    @show i;
    SPAC.spac!(config, spac, FT(1));
end;


gs_sl = spac.plant.leaves[3].flux.state.g_H₂O_s[1:end-1];
gs_sh = spac.plant.leaves[3].flux.state.g_H₂O_s[end];
@show maximum(gs_sl);
@show minimum(gs_sl);
@show gs_sh;

an_sl = spac.plant.leaves[3].flux.auxil.a_n[1:end-1];
an_sh = spac.plant.leaves[3].flux.auxil.a_n[end];
@show maximum(an_sl);
@show minimum(an_sl);
@show an_sh;

et_sl, et_sh = Emerald.Land.ET_LEAF(config, spac);
et_sl_end = et_sl[:, :, end-2];
et_sh_end = et_sh[end-2];
@show maximum(et_sl_end);
@show minimum(et_sl_end);
@show et_sh_end;

wue_sl = an_sl ./ et_sl_end[:];
wue_sh = an_sh ./ et_sh_end;
@show maximum(wue_sl);
@show minimum(wue_sl);
@show wue_sh;

iwue_sl = an_sl ./ gs_sl;
iwue_sh = an_sh ./ gs_sh;
@show maximum(iwue_sl);
@show minimum(iwue_sl);
@show iwue_sh;


using DataFrames
df = DataFrame();
df[!,:PPAR] = spac.plant.leaves[3].flux.auxil.ppar;
df[!,:CI] = spac.plant.leaves[3].flux.auxil.p_CO₂_i;
df[!,:E] = [et_sl_end[:]...;et_sh_end];
df[!,:AN] = spac.plant.leaves[3].flux.auxil.a_n;
df[!,:GS] = spac.plant.leaves[3].flux.state.g_H₂O_s;
save_csv!(df, "test.csv")
