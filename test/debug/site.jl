import Emerald.Land as ELAND


settings = ELAND.land_model_settings(mode = "default");
settings["C3_MODEL"] = "FvCB";
settings["MESSAGE_LEVEL"] = 2;
nt = ELAND.simulation!(settings, -10.5, -70.5, 2019);
