import Emerald.Land as ELAND


settings = ELAND.land_model_settings(mode = "default");
settings["C3_MODEL"] = "Jmax";
settings["MESSAGE_LEVEL"] = 2;
nt = ELAND.simulation!(settings, -10.5, -70.5, 2019);

settings = ELAND.land_model_settings(mode = "ash");
settings["C3_MODEL"] = "Jmax";
settings["MESSAGE_LEVEL"] = 2;
nt = ELAND.simulation!(settings, -10.5, -70.5, 2023);
