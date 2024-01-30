EUNISpath <- "/data/GEO/emodnethabitats/"

fn=paste0(EUNISpath,"EUSM_Arctic_Atlantic.R")
load(file=fn)

r=raster::raster(paste0(EUNISpath,"EUSM_Arctic_Atlantic.tif"))
mapview(r)
r <- crop(r, c(-20, 11, 40, 65))
mapview(r)

head(EUSMdf)   ## code gives the current label of all classes
unique(EUSMdf$Substrate)

substr_lvl <- data.frame(levels = unique(EUSMdf$Substrate),
                         code_substr = c(1:12))

energy_lvl <- data.frame(levels = c("High energy", "Moderate energy", "Low energy", "No energy information"),
                         code_ene = c(1:4))

nrow(EUSMdf)

EUSMdf <- EUSMdf %>% left_join(substr_lvl, by = join_by(Substrate == levels))
EUSMdf <- EUSMdf %>% left_join(energy_lvl, by = join_by(Energy == levels))
head(EUSMdf)

r_substr <- reclassify(r, cbind(EUSMdf$code, EUSMdf$code_substr))
r_ene <- reclassify(r, cbind(EUSMdf$code, EUSMdf$code_ene))

mapview(r_ene)
mapview(r_substr)

writeRaster(r_ene, "seabed_energy_EMODnet.tif", overwrite = TRUE)
writeRaster(r_substr, "seabed_substrate_EMODnet.tif", overwrite = TRUE)
