library(raster)
library(tidyverse)
library(terrainr)
library(sf)
library(sp)
library(stars)

setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/")

wf_poly <- sf::st_read("GIS/Herring/abiotic layers/Windmills EU - EMODnet/EMODnet_HA_Energy_WindFarms_20230605/EMODnet_HA_Energy_WindFarms_pg_20230605.shp",
                    layer = "EMODnet_HA_Energy_WindFarms_pg_20230605")
ex_r <- stack("Scripts data-driven approach/Herring/Input-output_files/2.Combined files/chlorophyll_concentration_chl.grd")[[1]]
unique(wf_poly$STATUS)

print(as.data.frame(wf_poly) %>%
  dplyr::group_by(STATUS) %>% 
  dplyr::distinct(YEAR) %>%
    arrange(STATUS), n = 50)

# keep only production

wf_poly2 <- wf_poly %>% 
  filter(STATUS == "Production",
         YEAR != 0)

st <- stack()
for (y in 2000:2020) {
  for (m in 1:12) {
    poly <- wf_poly2 %>% 
      filter(YEAR <= y)
    st_crs(poly) <- st_crs("+init=epsg:4326")
    poly_b <- poly %>%
      st_transform("+init=epsg:32631") %>%
      st_buffer(200) %>%
      st_transform("+init=epsg:4326")
    r <- rasterize(poly_b, ex_r, 1, background = 0)
    r <- mask(r, ex_r)
    r <- crop(r, ex_r)
    st <- stack(st, r)
  }
  print(y)
}

t_index <- tibble(year = rep(c(2000:2020), each = 12),
                  month = rep(seq(1:12),21))
names(st) <- paste(t_index$year, t_index$month, sep = "_")
setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Scripts data-driven approach/Herring/Input-output_files/0.Preprocessed files")
writeRaster(st, "windfarms.grd", overwrite = TRUE)
