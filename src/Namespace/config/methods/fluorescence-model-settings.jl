# van der Tol et al. (2013) Models of fluorescence and photosynthesis for interpreting measurements of solar-induced chlorophyll fluorescence
KNFluorescenceModelAll(FT) = KNFluorescenceModel{FT}(K_0 = 2.48, K_A = 2.83, K_B = 0.114)
KNFluorescenceModelDrought(FT) = KNFluorescenceModel{FT}(K_0 = 5.01, K_A = 1.93, K_B = 10);


# Han et al. (2022) The physiological basis for estimating photosynthesis from Chla fluorescence
QLFluorescenceModelC3(FT) = QLFluorescenceModel{FT}(K_B = 0.95e-3 / 0.85);
QLFluorescenceModelC4(FT) = QLFluorescenceModel{FT}(K_B = 0.63e-3 / 0.85);
QLFluorescenceModelHanC3(FT) = QLFluorescenceModelHan{FT}(K_A = 0.8, K_B = 0.95e-3 / 0.85);
QLFluorescenceModelHanC4(FT) = QLFluorescenceModelHan{FT}(K_A = 0.83, K_B = 0.63e-3 / 0.85);
