
# Sources
#     Lavigne and Ryan (1997) Growth and maintenance respiration rates of aspen, blackspruce and jack pine stems at northern and southern BOREAS sites
#     Bernacchi et al. (2001) Improved temperature response functions for models of Rubisco‐limited photosynthesis
#     Boyd et al. (2001) Temperature responses of C4 photosynthesis: biochemical analysis of Rubisco, phosphoenolpyruvate carboxylase, and carbonic anhydrase in Setaria viridis
#     Leuning (2002) Temperature dependence of two parameters in a photosynthesis model
#     Kattge et al. (2007) Temperature acclimation in a biochemical model of photosynthesis: a reanalysis of data from 36 species
#     Sperry et al. (2019) The impact of rising CO2 and acclimation on the response of US forests to global warming
#     Johnson et al. (2021) The limiting factors and regulatory processes that control the environmental responses of C3, C3–C4 intermediate, and C4 photosynthesis
#     CLM5 Documentation. Chapter 9 Page 106
KcTDBernacchi(FT)          = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 41.0264925, ΔHA = 79430.0);
KcTDCLM(FT)                = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 40.49     , ΔHA = 79430.0);
KoTDBernacchi(FT)          = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 28208.88  , ΔHA = 36380.0);
KoTDCLM(FT)                = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 27840.0   , ΔHA = 36380.0);
KpepTDBoyd(FT)             = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 16.0      , ΔHA = 36300.0);
KqTDJohnson(FT)            = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 300       , ΔHA = 37000.0);
RespirationTDBernacchi(FT) = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = NaN       , ΔHA = 46390.0);
VcmaxTDBernacchi(FT)       = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = NaN       , ΔHA = 65330.0);
VomaxTDBernacchi(FT)       = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = NaN       , ΔHA = 60110.0);
ΓStarTDBernacchi(FT)       = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 4.33164375, ΔHA = 37830.0);
ΓStarTDCLM(FT)             = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 4.275     , ΔHA = 37830.0);

JmaxTDBernacchi(FT)                 = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 57500.0, ΔHD = 439000.0, ΔSV = 1400.0);
JmaxTDCLM(FT, t::Number = T₂₅())    = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 50000.0, ΔHD = 200000.0, ΔSV = 659.70 - 0.75 * (t - T₀(FT)) );
JmaxTDLeuning(FT)                   = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 50300.0, ΔHD = 152044.0, ΔSV = 495.0 );
RespirationTDCLMC3(FT)              = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 46390.0, ΔHD = 150650.0, ΔSV = 490.0 );
VcmaxTDCLMC3(FT, t::Number = T₂₅()) = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 72000.0, ΔHD = 200000.0, ΔSV = 668.39 - 1.07 * (t - T₀(FT)) );
VcmaxTDLeuning(FT)                  = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 73637.0, ΔHD = 149252.0, ΔSV = 486.0 );
VpmaxTDBoyd(FT)                     = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN , ΔHA = 94800.0, ΔHD = 73300.0 , ΔSV = 250.0 );
ηCTDJohnson(FT)                     = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 1.0 , ΔHA = 0.0    , ΔHD = 220000.0, ΔSV = 710.0 );
ηLTDJohnson(FT)                     = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 0.75, ΔHA = 0.0    , ΔHD = 220000.0, ΔSV = 710.0 );
ηCTDWang(FT)                        = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 1.0 , ΔHA = 0.0    , ΔHD = 225100.0, ΔSV = 710.0 );
ηLTDWang(FT)                        = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 0.75, ΔHA = 0.0    , ΔHD = 225100.0, ΔSV = 710.0 );

Q10TDAngiosperm(FT) = Q10{FT}(Q_10 = 1.4, T_REF = T₂₅(FT), VAL_REF = 2 * 0.0140 / 8760 * 1000); # μmol CO2 mol⁻¹ C biomass s⁻¹
Q10TDGymnosperm(FT) = Q10{FT}(Q_10 = 1.7, T_REF = T₂₅(FT), VAL_REF = 2 * 0.0425 / 8760 * 1000); # μmol CO2 mol⁻¹ C biomass s⁻¹
Q10TDKpepCLM(FT)    = Q10{FT}(Q_10 = 2.0, T_REF = T₂₅(FT), VAL_REF = 0.2);

RespirationTDCLMC4(FT) = Q10PeakHT{FT}(Q_10 = 2.0, T_REF = T₂₅(FT), VAL_REF = NaN, ΔT_REF = 328.15, ΔT_SLOPE = 1.3);

VcmaxTDCLMC4(FT) = Q10PeakLTHT{FT}(Q_10 = 2.0, T_REF = T₂₅(FT), VAL_REF = NaN, ΔHT_REF = 313.15, ΔHT_SLOPE = 0.3, ΔLT_REF = 288.15, ΔLT_SLOPE = 0.2);


# New parameters for the temperature dependency based on fitting A-Ci curves I collected
# TODO: make it default in the future after the paper is accepted
ΓStarTDWang2024(FT) = Arrhenius{FT}(T_REF = T₂₅(FT), VAL_REF = 4.56, ΔHA = 11800.0);

JmaxTDWang2024(FT, t::Number = T₂₅())  = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN   , ΔHA = 50000, ΔHD = 201000, ΔSV = 659.70 - 0.75 * (t - T₀(FT)));
KqTDWang2024(FT)                       = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 300   , ΔHA = 21900, ΔHD = 232000, ΔSV = 700);
VcmaxTDWang2024(FT, t::Number = T₂₅()) = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = NaN   , ΔHA = 63000, ΔHD = 204000, ΔSV = 668.39 - 1.07 * (t - T₀(FT)));
ηCTDWang2024(FT)                       = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 2*3/14, ΔHA = 21900, ΔHD = 232000, ΔSV = 700);
ηLTDWang2024(FT)                       = ArrheniusPeak{FT}(T_REF = T₂₅(FT), VAL_REF = 3*3/14, ΔHA = 21900, ΔHD = 232000, ΔSV = 700);
