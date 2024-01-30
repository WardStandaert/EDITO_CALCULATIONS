library(tidyverse)
setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling")
haul_data <- read.csv("GIS/Herring/biotic layers/Herring Datras/ICESDataPortalDownload_DATRAS_93345b25-d93a-47fd-bd57-2915b03ce7f6/HL_0703271736.csv")
location <- read.csv("GIS/Herring/biotic layers/Herring Datras/ICESDataPortalDownload_DATRAS_93345b25-d93a-47fd-bd57-2915b03ce7f6/location_v2.csv")
glimpse(haul_data)
glimpse(location)

world <- map_data("world")
europe <- subset(world, region %in% c("Albania", "Andorra", "Armenia", "Austria", "Azerbaijan",
                                      "Belarus", "Belgium", "Bosnia and Herzegovina", "Bulgaria",
                                      "Croatia", "Cyprus", "Czechia","Denmark","Estonia","Finland",
                                      "France","Georgia", "Germany", "Greece","Hungary","Iceland", 
                                      "Ireland", "Italy","Kazakhstan", "Kosovo", "Latvia","Liechtenstein", 
                                      "Lithuania", "Luxembourg","Malta","Moldova","Monaco","Montenegro",
                                      "Macedonia", "Netherlands","Norway","Poland","Portugal","Romania",
                                      "Russia","San Marino","Serbia","Slovakia","Slovenia","Spain",
                                      "Sweden","Switzerland","Turkey","Ukraine","UK","Vatican"))

#### FISH DATA PROCESSING ####
#filter species based on their AphiaID and year of sampling
haul_data2 <- haul_data %>%
  filter(ValidAphiaID == 126417 |
           ValidAphiaID == 126413 |
           ValidAphiaID == 127023 |
           ValidAphiaID == 126975,
         Year >= 2000 & Year <= 2020)
unique(haul_data2$ScientificName_WoRMS)

#explore for some measure of spawning
haul_data2 %>% group_by(ScientificName_WoRMS, DevStage) %>%
  summarize(n = n())
haul_data2 %>% group_by(ScientificName_WoRMS, SubWgt) %>%
  summarize(n = n()) %>%
  arrange(desc(n))
haul_data2 %>% group_by(ScientificName_WoRMS, LngtClass) %>%
  summarize(n = n()) %>%
  arrange(desc(n))
#no dev stages

haul_data2 %>% filter(!(LngtClass == -9 | LngtClass == -9)) %>%
  group_by(ScientificName_WoRMS) %>%
  count

#remove observations without length measurement
hauls_w_len <- haul_data2 %>% filter(!(LngtClass == -9 | LngtClass == -9))








#### HAUL LOCATION PROCESSING ####
# convert to numeric + add columns to indicate validity of lon-lat pairs
location2 <- location %>%
  mutate(HaulLong = as.numeric(HaulLong),
         HaulLat = as.numeric(HaulLat),
         ShootLong = as.numeric(ShootLong),
         ShootLat = as.numeric(ShootLat),
         Lon = rep(NA, nrow(location)),
         Lat = rep(NA, nrow(location)))

# there are two values representing NA: NA and -9
# convert all NA data to -9 so that there is only one value
location2$HaulLong[which(is.na(location2$HaulLong))] <- -9
location2$HaulLat[which(is.na(location2$HaulLat))] <- -9
location2$ShootLong[which(is.na(location2$ShootLong))] <- -9
location2$ShootLat[which(is.na(location2$ShootLat))] <- -9

which(colnames(location2) %in% c("ShootLong", "ShootLat", "HaulLong", "HaulLat"))

# retrieve valid lon and lat values
valid_coords <- function(vector) {
  # 1. shootcoords
  ShootCoord <- as.numeric(c(vector[20], vector[19]))
  ShootVal <- ifelse(-9 %in% ShootCoord, "not_valid", "valid")
  
  # 2. haulcoords
  HaulCoord <- as.numeric(c(vector[22], vector[21]))
  HaulVal <- ifelse(-9 %in% HaulCoord, "not_valid", "valid")
  
  # create output
  if(ShootVal == "not_valid" & HaulVal == "not_valid") {
    Lon <- NA
    Lat <- NA
  } else if (ShootVal == "not_valid" & HaulVal == "valid") {
    Lon <- HaulCoord[1]
    Lat <- HaulCoord[2]
  } else if (ShootVal == "valid" & HaulVal == "not_valid") {
    Lon <- ShootCoord[1]
    Lat <- ShootCoord[2]
  } else {
    Lon <- mean(ShootCoord[1], HaulCoord[1])
    Lat <- mean(ShootCoord[2], HaulCoord[2])
  }
  return(as.numeric(cbind(Lon, Lat)))
}

LonLat <- apply(location2, 1, valid_coords)
location2$Lon <- t(LonLat)[,1]
location2$Lat <- t(LonLat)[,2]

# combine filter on location
location3 <- location2 %>%
  filter(Lon > -20 & Lon < 30,
         Lat > 40 & Lat < 65)

location3 %>% group_by(Lon) %>%
  summarize(n = n()) %>% 
  arrange(desc(n))
location3 %>% group_by(Lat) %>%
  summarize(n = n()) %>% 
  arrange(desc(n))

#remove lon's that equal 1 (seems like error)
location4 <- location3 %>%
  filter(Lon != 1)

location4 %>% 
  group_by(StNo) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

stations_loc_summary <- location4 %>% 
  filter(StNo != -9) %>%
  group_by(StNo, Year) %>%
  summarize(n = n(),
            Lon_min = min(Lon),
            Lon_mean = mean(Lon),
            Lon_max = max(Lon),
            Lon_range = Lon_min - Lon_max,
            Lat_min = min(Lon),
            Lat_mean = mean(Lon),
            Lat_max = max(Lon),
            Lat_range = Lat_min - Lat_max)

stations_loc_summary %>% group_by(Lon_range, Lat_range) %>%
  summarize(n = n()) %>%
  arrange(desc(Lon_range))

stations_w_consistent_location <- location4 %>% 
  filter(StNo != -9) %>%
  group_by(StNo, Year) %>%
  summarize(n = n(),
            Lon_min = min(Lon),
            Lon_mean = mean(Lon),
            Lon_max = max(Lon),
            Lon_range = Lon_min - Lon_max,
            Lat_min = min(Lat),
            Lat_mean = mean(Lat),
            Lat_max = max(Lat),
            Lat_range = Lat_min - Lat_max) %>%
  filter(abs(Lon_range) <= 30 & abs(Lat_range) <= 30)


ggplot(location4, aes(Lon, Lat)) +
  geom_point() + 
  geom_polygon(data = europe, aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  theme_void() +
  coord_fixed(ratio=1.6, xlim = c(-10, 30),
              ylim = c(40,65)) 


#### MERGING FISH AND LOCATION DATA ####

df <- hauls_w_len %>% filter(trim(StNo) %in% stations_w_consistent_location$StNo)

hauls_w_len %>% filter(!(StNo %in% unique(df$StNo))) %>%
  group_by(StNo) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

df %>% group_by(ScientificName_WoRMS) %>%
  summarize(n = n())

hauls_w_len <- hauls_w_len %>%
  mutate(StNo = trim(StNo))

df <- left_join(hauls_w_len, select(stations_w_consistent_location, StNo, Lon = Lon_mean, Lat = Lat_mean), by = "StNo", multiple = "first")
df <- df %>% drop_na(Lon, Lat)

##remove duplicates
df_small <- df %>%
  filter(LngtClass <= 10)
df_big <- df %>%
  filter(LngtClass > 10)

df_small_u <- unique(data.frame(ScientificName_WoRMS = df_small$ScientificName_WoRMS, 
                                Lon = df_small$Lon, 
                                Lat = df_small$Lat))
df_big_u <- unique(data.frame(ScientificName_WoRMS = df_big$ScientificName_WoRMS, 
                              Lon = df_big$Lon, 
                              Lat = df_big$Lat))


ggplot(df_big_u, aes(Lon, Lat, color = ScientificName_WoRMS)) +
  geom_point() + 
  geom_polygon(data = europe, aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  theme_void() +
  coord_fixed(ratio=1.6, xlim = c(-10, 30),
              ylim = c(40,65)) 

ggplot(df_small_u, aes(Lon, Lat, color = ScientificName_WoRMS)) +
  geom_point() + 
  geom_polygon(data = europe, aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  theme_void() +
  coord_fixed(ratio=1.6, xlim = c(-10, 30),
              ylim = c(40,65)) 

write.csv(df_big_u, "Scripts data-driven approach/Herring/Input-output_files/2.Combined files/dummy_data_adult.csv")
write.csv(df_small_u, "Scripts data-driven approach/Herring/Input-output_files/2.Combined files/dummy_data_larvae.csv")
