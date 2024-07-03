# 1. Prepare data extraction ----
setwd("/home/onyxia/work/EDITO_CALCULATIONS/demonstrators/Ward/zarr_extraction/")

# load functions, but also cached stac catalog named EDITOSTAC
source("editoTools.R")
options("outputdebug"=c('L','M'))

# the file to process
# datafile = "./data-raw/herring_test_set.csv"
datafile = "./herring_zarr_extraction_test_output.csv"
points = read.delim(datafile, sep=",") %>% sample_n(100)
glimpse(points)

points <- points %>%
  dplyr::select(Longitude, Latitude, Year, Month, Time, SST, SSS, windfarms, ZooPl, Phyto, seabed_energy, seabed_substrate, depth)

#the requested timestep resolution of the dataset in milliseconds
#in this case we work with monthly data (1 month = 30.436875*24*3600*1000 = 2629746000 milliseconds)
timeSteps=c(2629746000)


#the requested parameters, names in the stac catalogue table field 'par' .. see unique(stacCatalogue$par) for a list 
parameters = list("zooc"= c("par" = "zooc",
                            "fun" = "mean",
                            "buffer" = "10000",
                            "convert_from_timestep" = 86400000),
                  "phyc"= c("par" = "phyc",
                            "fun" = "mean",
                            "buffer" = "10000"),
                  "thetao"= c("par" = "thetao",
                              "fun" = "mean",
                              "buffer" = "10000"),
                  "so"= c("par" = "so",
                          "fun" = "mean",
                          "buffer" = "10000"))

# parameters = list("Substrate"= c("par" = "Substrate",
#                                  "fun" = "table",
#                                  "buffer" = "10000"))

# parameters= list("Energy"= c("par" = "Energy", 
#                              "fun" = "most_freq", 
#                              "buffer" = "10000"),
#                  "Substrate"= c("par" = "Substrate",
#                                 "fun" = "most_freq",
#                                 "buffer"= "10000"))


#check if they all exist
for ( parameter in parameters) {
  param = ifelse(length(parameter) == 1, parameter, parameter["par"])
  if(! param %in% unique(EDITOSTAC$par) ) 
    { dbl("Unknown parameter ", param)
  }
}

# 2. Extract data ----
#add verbose= anything to get additional info on the positions ( par_x, par_y, par_z ) and time (par_t) found in the zarr files
source("editoTools.R")
enhanced_DF = enhanceDF(inputPoints = points,
                        requestedParameters = parameters, 
                        requestedTimeSteps = NA, 
                        stacCatalogue = EDITOSTAC, 
                        verbose="",
                        select_layers = rep(1,length(parameters)),
                        atDepth = 20)


#try to add new zarr file to catalog manually (downsampled elevation)
#start by loading the editoTools, this also loads the EDITOSTAC catalog
source("editoTools.R")
glimpse(EDITOSTAC)
#add the new bathymetry file to the catalog
url <- "https://minio.lab.dive.edito.eu/oidc-willemboone/bathymetry_ward.zarr"
EDITOSTAC_ed <- EDITOSTAC %>% add_row(par = "elevation",
                      href = url,
                      latmin = 48,
                      latmax = 62,
                      lonmin = -12,
                      lonmax = 10,
                      start_datetime = "1900-01-01",
                      end_datetime = "3000-01-01",
                      title = "downsampled_bathymetry",
                      categories = 0,
                      catalogue = "EMODNET",
                      chunktype = "chunked")

parameters = list("elevation"= c("par" = "elevation",
                                 "fun" = "mean",
                                 "buffer" = "18000"))

enhanced_DF_depth = enhanceDF(inputPoints = points,
                        requestedParameters = parameters, 
                        requestedTimeSteps = NA, 
                        stacCatalogue = EDITOSTAC_ed, 
                        verbose="on",
                        atDepth = 20)


#1,1,1,1
#always opt for timeChunked
glimpse(enhanced_DF_depth)

## trials on extracting seabed habitat data
# edit EDITOSTAC twice (should be solved)
EDITOSTAC[which(EDITOSTAC$par == "EUNIS2019C" & EDITOSTAC$asset == "Zarr"),c("latmin","latmax","lonmin","lonmax")] <- matrix(data = c(48,62,-12,10), nrow = 1)
EDITOSTAC[which(EDITOSTAC$par == "EUNIS2019C" & EDITOSTAC$asset == "Zarr"),"par"] <- "eunis_seabed_habitat_class_2019"
r <- getRasterSlice(requestedParameter = "eunis_seabed_habitat_class_2019", lon_min = -12, lon_max = 10, lat_min = 48, lat_max = 62, stacCatalogue = EDITOSTAC)
plot(r)



#test on selecting specific categorical variable
parameters = list("seabed_energy"= c("par" = "seabed_energy",
                              "fun" = "table",
                              "buffer" = "10000"))
EDITOSTAC[which(EDITOSTAC$par == "Energy" & EDITOSTAC$asset == "Zarr"),c("latmin","latmax","lonmin","lonmax")] <- matrix(data = c(48,62,-12,10), nrow = 1)
EDITOSTAC[which(EDITOSTAC$par == "Energy" & EDITOSTAC$asset == "Zarr"),c("latmin","latmax","lonmin","lonmax")] <- matrix(data = c(48,62,-12,10), nrow = 1)
EDITOSTAC[which(EDITOSTAC$par == "Energy" & EDITOSTAC$asset == "Zarr"),"par"] <- "seabed_energy"
enhanced_DF_depth = enhanceDF(inputPoints = points,
                              requestedParameters = parameters, 
                              requestedTimeSteps = NA, 
                              stacCatalogue = EDITOSTAC, 
                              verbose="on",
                              atDepth = 20)

r <- getRasterSlice(requestedParameter = "seabed_energy",
                    lon_min = -12,
                    lon_max = 10,
                    lat_min = 48,
                    lat_max = 62,
                    requestedTimeSteps = NA,
                    date = "2020-01-01",
                    stacCatalogue = EDITOSTAC)
plot(r)

#does not work if no period is provided. would work once monthly data is not represented as NA anymore
new_pts <- tibble(Longitude = c(0,2,4,6,8),
                  Latitude = c(50,50,50,50,50),
                  Time = rep(as.POSIXct("2020-01-01",tz = "UTC")))

enhanced_DF = enhanceDF(inputPoints = new_pts,
                              requestedParameters = list("seabed_energy"= c("par" = "seabed_energy",
                                                                            "fun" = "most_frequent",
                                                                            "buffer" = "50000")), 
                              requestedTimeSteps = NA, 
                              stacCatalogue = EDITOSTAC, 
                              verbose="on",
                              atDepth = 0)
enhanced_DF$seabed_energy
# [1] 3 3 4 4 2

#try an exact match with category 3 (e.g. species likes to be close to corals)
enhanced_DF = enhanceDF(inputPoints = new_pts,
                        requestedParameters = list("seabed_energy"= c("par" = "seabed_energy",
                                                                      "fun" = "exact",
                                                                      "category" = "3",
                                                                      "buffer" = "50000")), 
                        requestedTimeSteps = NA, 
                        stacCatalogue = EDITOSTAC, 
                        verbose="on",
                        atDepth = 0)
enhanced_DF$seabed_energy
# [1] 3 3 4 4 2

## Compare with original extraction  ----
### Numerical variables ----
par(mfrow = c(3,2))
plot(enhanced_DF_depth$SST, enhanced_DF_depth$thetao, cex.lab=1.3,
     xlim = c(min(enhanced_DF_depth$SST, enhanced_DF_depth$thetao, na.rm = T), max(enhanced_DF_depth$SST, enhanced_DF_depth$thetao, na.rm = T)),
     ylim = c(min(enhanced_DF_depth$SST, enhanced_DF_depth$thetao, na.rm = T), max(enhanced_DF_depth$SST, enhanced_DF_depth$thetao, na.rm = T)))
plot(enhanced_DF_depth$SSS, enhanced_DF_depth$so, cex.lab=1.3,
     xlim = c(min(enhanced_DF_depth$SSS, enhanced_DF_depth$so, na.rm = T), max(enhanced_DF_depth$SSS, enhanced_DF_depth$so, na.rm = T)),
     ylim = c(min(enhanced_DF_depth$SSS, enhanced_DF_depth$so, na.rm = T), max(enhanced_DF_depth$SSS, enhanced_DF_depth$so, na.rm = T)))
plot(enhanced_DF_depth$Phyto, enhanced_DF_depth$phyc, cex.lab=1.3,
     xlim = c(min(enhanced_DF_depth$Phyto, enhanced_DF_depth$phyc, na.rm = T), max(enhanced_DF_depth$Phyto, enhanced_DF_depth$phyc, na.rm = T)),
     ylim = c(min(enhanced_DF_depth$Phyto, enhanced_DF_depth$phyc, na.rm = T), max(enhanced_DF_depth$Phyto, enhanced_DF_depth$phyc, na.rm = T)))
plot(enhanced_DF_depth$ZooPl, enhanced_DF_depth$zooc, cex.lab=1.3,
     xlim = c(min(enhanced_DF_depth$ZooPl, enhanced_DF_depth$zooc, na.rm = T), max(enhanced_DF_depth$ZooPl, enhanced_DF_depth$zooc, na.rm = T)),
     ylim = c(min(enhanced_DF_depth$ZooPl, enhanced_DF_depth$zooc, na.rm = T), max(enhanced_DF_depth$ZooPl, enhanced_DF_depth$zooc, na.rm = T)))
plot(enhanced_DF_depth$depth, enhanced_DF_depth$elevation, cex.lab=1.3,
     xlim = c(min(enhanced_DF_depth$depth, enhanced_DF_depth$elevation, na.rm = T), max(enhanced_DF_depth$depth, enhanced_DF_depth$elevation, na.rm = T)),
     ylim = c(min(enhanced_DF_depth$depth, enhanced_DF_depth$elevation, na.rm = T), max(enhanced_DF_depth$depth, enhanced_DF_depth$elevation, na.rm = T)))


# 3. Check extract outcome ----
## Check difference in time  ----

difftime(enhanced_DF_depth$Time, enhanced_DF_depth$thetao_t, units = "days")
difftime(enhanced_DF_depth$Time, enhanced_DF_depth$so_t, units = "days")
difftime(enhanced_DF_depth$Time, enhanced_DF_depth$zooc_t, units = "days")
difftime(enhanced_DF_depth$Time, enhanced_DF_depth$phyc_t, units = "days")

## Check difference in space  ----
library(geosphere)


enhanced_DF_depth$thetao_dist <- apply(enhanced_DF_depth, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  thetao_x <- as.numeric(row["thetao_x"])
  thetao_y <- as.numeric(row["thetao_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(thetao_x, thetao_y), fun = distHaversine)
})

enhanced_DF_depth$so_dist <- apply(enhanced_DF_depth, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  so_x <- as.numeric(row["so_x"])
  so_y <- as.numeric(row["so_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(so_x, so_y), fun = distHaversine)
})

enhanced_DF_depth$zooc_dist <- apply(enhanced_DF_depth, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  zooc_x <- as.numeric(row["zooc_x"])
  zooc_y <- as.numeric(row["zooc_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(zooc_x, zooc_y), fun = distHaversine)
})

enhanced_DF_depth$phyc_dist <- apply(enhanced_DF_depth, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  phyc_x <- as.numeric(row["phyc_x"])
  phyc_y <- as.numeric(row["phyc_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(phyc_x, phyc_y), fun = distHaversine)
})

enhanced_DF_depth$elevation_dist <- apply(enhanced_DF_depth, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  elevation_x <- as.numeric(row["elevation_x"])
  elevation_y <- as.numeric(row["elevation_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(elevation_x, elevation_y), fun = distHaversine)
})

enhanced_DF_depth$Substrate_dist <- apply(enhanced_DF_depth, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  Substrate_x <- as.numeric(row["Substrate_x"])
  Substrate_y <- as.numeric(row["Substrate_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(Substrate_x, Substrate_y), fun = distHaversine)
})

enhanced_DF_depth$Energy_dist <- apply(enhanced_DF_depth, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  Energy_x <- as.numeric(row["Energy_x"])
  Energy_y <- as.numeric(row["Energy_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(Energy_x, Energy_y), fun = distHaversine)
})

glimpse(enhanced_DF_depth)

par(mfrow = c(3,3))
hist(enhanced_DF_depth$thetao_dist)
hist(enhanced_DF_depth$so_dist)
hist(enhanced_DF_depth$phyc_dist)
hist(enhanced_DF_depth$zooc_dist)
hist(enhanced_DF_depth$elevation_dist)
hist(enhanced_DF_depth$Substrate_dist)
hist(enhanced_DF_depth$Energy_dist)


### Categorical variables ----
substr_lvl <- tibble(sub_char = c("Fine mud", "Sand", "Muddy sand", "Mixed sediment",
                                  "Coarse substrate","Sandy mud or Muddy sand", "Seabed",
                                  "Rock or other hard substrata","Sandy mud", "Sandy mud or Muddy sand ",
                                  "Sediment","Fine mud or Sandy mud or Muddy sand"),
                     seabed_substrate = c(1:12))
energy_lvl <- tibble(ene_char = c("High energy", "Moderate energy", "Low energy", "No energy information"),
                     seabed_energy = c(1:4))

enhanced_DF_depth <- enhanced_DF_depth %>%
  left_join(substr_lvl, by = c("seabed_substrate")) %>%
  left_join(energy_lvl, by = c("seabed_energy"))

glimpse(enhanced_DF_depth)

sum(enhanced_DF_depth$Energy_Description == enhanced_DF_depth$ene_char) / length(enhanced_DF_depth$Energy_Description)
# 94% match
sum(enhanced_DF_depth$Substrate_Description == enhanced_DF_depth$sub_char) / length(enhanced_DF_depth$Substrate_Description)
# 82% match


# 4. Extract raster slice from .zarr ----
source("editoTools.R")
options("outputdebug"=c('L','M'))
load(file = "editostacv2.par")

#the requested timestep resolution of the dataset in milliseconds
#in this case we work with monthly data (1 month = 30.436875*24*3600*1000 = 2629746000 milliseconds)
timeSteps=c(2629746000)

EDITOSTAC[which(EDITOSTAC$par == "elevation" & EDITOSTAC$asset == "Zarr"),c("latmin","latmax","lonmin","lonmax")] <- matrix(data = c(25,85,-36,43), nrow = 1)
r <- getRasterSlice(requestedParameter = "elevation",
                    lon_min = -10,
                    lon_max = 10,
                    lat_min = 50,
                    lat_max = 55,
                    requestedTimeSteps = NA,
                    date = "2020-01-01",
                    stacCatalogue = EDITOSTAC)

# with predefined user input
r <- getRasterSlice(requestedParameter = "thetao",
                    lon_min = -13,
                    lon_max = 10,
                    lat_min = 40,
                    lat_max = 60,
                    requestedTimeSteps = NA,
                    date = "2020-01-01",
                    stacCatalogue = EDITOSTAC,
                    select_layers = 1)

# with predefined user input
EDITOSTAC[which(EDITOSTAC$par == "Energy" & EDITOSTAC$asset == "Zarr"),c("latmin","latmax","lonmin","lonmax")] <- matrix(data = c(20,80,-40,40), nrow = 1)
EDITOSTAC[which(EDITOSTAC$par == "Energy" & EDITOSTAC$asset == "Zarr"),"par"] <- "seabed_energy"
r <- getRasterSlice(requestedParameter = "seabed_energy",
                    lon_min = -13,
                    lon_max = 10,
                    lat_min = 40,
                    lat_max = 60,
                    requestedTimeSteps = NA,
                    stacCatalogue = EDITOSTAC)

plot(r)

#does not work
extract(r, terra::buffer(vect(cbind(c(0,0),c(40,50)), crs="+proj=longlat"), mean, width = 100))
#works
tab = table(terra::extract(r, terra::buffer(vect(cbind(c(0,0),c(40,50)), crs="+proj=longlat"), width = 10000)))
round(tab / rowSums(tab), 2)
