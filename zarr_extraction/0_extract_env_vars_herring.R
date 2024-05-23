# 1. Prepare data extraction ----
setwd("/home/onyxia/work/EDITO_CALCULATIONS/zarr_extraction/")

# load functions, but also cached stac catalog named stacCatalog
source("editoTools.R")
options("outputdebug"=c('L','M'))

#the cached stacCatalog is called 'EDITOSTAC'
load(file = "./data-raw/editostacv2.par")

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
parameters = list("thetao"= c("par" = "thetao",
                              "fun" = "mean",
                              "buffer" = "10000"),
                  "so"= c("par" = "so",
                          "fun" = "mean",
                          "buffer" = "10000"),
                  "zooc"= c("par" = "zooc",
                            "fun" = "mean",
                            "buffer" = "10000",
                            "convert_from_timestep" = 86400000),
                  "phyc"= c("par" = "phyc",
                            "fun" = "mean",
                            "buffer" = "10000"))

parameters = list("zooc"= c("par" = "zooc",
                            "fun" = "mean",
                            "buffer" = "10000",
                            "convert_from_timestep" = 86400000),
                  "phyc"= c("par" = "phyc",
                            "fun" = "mean",
                            "buffer" = "10000"))

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
                        verbose="on",
                        select_layers = rep(1,length(parameters)),
                        atDepth = 20)


#1,1,1,1
#always opt for timeChunked
glimpse(enhanced_DF)


## Compare with original extraction  ----
### Numerical variables ----
par(mfrow = c(3,2))
plot(enhanced_DF$SST, enhanced_DF$thetao, cex.lab=1.3,
     xlim = c(min(enhanced_DF$SST, enhanced_DF$thetao, na.rm = T), max(enhanced_DF$SST, enhanced_DF$thetao, na.rm = T)),
     ylim = c(min(enhanced_DF$SST, enhanced_DF$thetao, na.rm = T), max(enhanced_DF$SST, enhanced_DF$thetao, na.rm = T)))
plot(enhanced_DF$SSS, enhanced_DF$so, cex.lab=1.3,
     xlim = c(min(enhanced_DF$SSS, enhanced_DF$so, na.rm = T), max(enhanced_DF$SSS, enhanced_DF$so, na.rm = T)),
     ylim = c(min(enhanced_DF$SSS, enhanced_DF$so, na.rm = T), max(enhanced_DF$SSS, enhanced_DF$so, na.rm = T)))
plot(enhanced_DF$Phyto, enhanced_DF$phyc, cex.lab=1.3,
     xlim = c(min(enhanced_DF$Phyto, enhanced_DF$phyc, na.rm = T), max(enhanced_DF$Phyto, enhanced_DF$phyc, na.rm = T)),
     ylim = c(min(enhanced_DF$Phyto, enhanced_DF$phyc, na.rm = T), max(enhanced_DF$Phyto, enhanced_DF$phyc, na.rm = T)))
plot(enhanced_DF$ZooPl, enhanced_DF$zooc, cex.lab=1.3,
     xlim = c(min(enhanced_DF$ZooPl, enhanced_DF$zooc, na.rm = T), max(enhanced_DF$ZooPl, enhanced_DF$zooc, na.rm = T)),
     ylim = c(min(enhanced_DF$ZooPl, enhanced_DF$zooc, na.rm = T), max(enhanced_DF$ZooPl, enhanced_DF$zooc, na.rm = T)))
# plot(enhanced_DF$depth, enhanced_DF$elevation, cex.lab=1.3,
#      xlim = c(min(enhanced_DF$depth, enhanced_DF$elevation, na.rm = T), max(enhanced_DF$depth, enhanced_DF$elevation, na.rm = T)),
#      ylim = c(min(enhanced_DF$depth, enhanced_DF$elevation, na.rm = T), max(enhanced_DF$depth, enhanced_DF$elevation, na.rm = T)))


# 3. Check extract outcome ----
## Check difference in time  ----

difftime(enhanced_DF$Time, enhanced_DF$thetao_t, units = "days")
difftime(enhanced_DF$Time, enhanced_DF$so_t, units = "days")
difftime(enhanced_DF$Time, enhanced_DF$zooc_t, units = "days")
difftime(enhanced_DF$Time, enhanced_DF$phyc_t, units = "days")

## Check difference in space  ----
library(geosphere)


enhanced_DF$thetao_dist <- apply(enhanced_DF, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  thetao_x <- as.numeric(row["thetao_x"])
  thetao_y <- as.numeric(row["thetao_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(thetao_x, thetao_y), fun = distHaversine)
})

enhanced_DF$so_dist <- apply(enhanced_DF, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  so_x <- as.numeric(row["so_x"])
  so_y <- as.numeric(row["so_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(so_x, so_y), fun = distHaversine)
})

enhanced_DF$zooc_dist <- apply(enhanced_DF, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  zooc_x <- as.numeric(row["zooc_x"])
  zooc_y <- as.numeric(row["zooc_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(zooc_x, zooc_y), fun = distHaversine)
})

enhanced_DF$phyc_dist <- apply(enhanced_DF, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  phyc_x <- as.numeric(row["phyc_x"])
  phyc_y <- as.numeric(row["phyc_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(phyc_x, phyc_y), fun = distHaversine)
})

enhanced_DF$elevation_dist <- apply(enhanced_DF, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  elevation_x <- as.numeric(row["elevation_x"])
  elevation_y <- as.numeric(row["elevation_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(elevation_x, elevation_y), fun = distHaversine)
})

enhanced_DF$Substrate_dist <- apply(enhanced_DF, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  Substrate_x <- as.numeric(row["Substrate_x"])
  Substrate_y <- as.numeric(row["Substrate_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(Substrate_x, Substrate_y), fun = distHaversine)
})

enhanced_DF$Energy_dist <- apply(enhanced_DF, 1, function(row) {
  # Convert elements to numeric
  Longitude <- as.numeric(row["Longitude"])
  Latitude <- as.numeric(row["Latitude"])
  Energy_x <- as.numeric(row["Energy_x"])
  Energy_y <- as.numeric(row["Energy_y"])
  
  # Calculate distance
  distm(c(Longitude, Latitude), c(Energy_x, Energy_y), fun = distHaversine)
})

glimpse(enhanced_DF)

par(mfrow = c(3,3))
hist(enhanced_DF$thetao_dist)
hist(enhanced_DF$so_dist)
hist(enhanced_DF$phyc_dist)
hist(enhanced_DF$zooc_dist)
# hist(enhanced_DF$elevation_dist)
hist(enhanced_DF$Substrate_dist)
hist(enhanced_DF$Energy_dist)


### Categorical variables ----
substr_lvl <- tibble(sub_char = c("Fine mud", "Sand", "Muddy sand", "Mixed sediment",
                                  "Coarse substrate","Sandy mud or Muddy sand", "Seabed",
                                  "Rock or other hard substrata","Sandy mud", "Sandy mud or Muddy sand ",
                                  "Sediment","Fine mud or Sandy mud or Muddy sand"),
                     seabed_substrate = c(1:12))
energy_lvl <- tibble(ene_char = c("High energy", "Moderate energy", "Low energy", "No energy information"),
                     seabed_energy = c(1:4))

enhanced_DF <- enhanced_DF %>%
  left_join(substr_lvl, by = c("seabed_substrate")) %>%
  left_join(energy_lvl, by = c("seabed_energy"))

glimpse(enhanced_DF)

sum(enhanced_DF$Energy_Description == enhanced_DF$ene_char) / length(enhanced_DF$Energy_Description)
# 94% match
sum(enhanced_DF$Substrate_Description == enhanced_DF$sub_char) / length(enhanced_DF$Substrate_Description)
# 82% match


# 4. Extract raster slice from .zarr ----
source("editoTools.R")
options("outputdebug"=c('L','M'))
load(file = "editostacv2.par")

#the requested timestep resolution of the dataset in milliseconds
#in this case we work with monthly data (1 month = 30.436875*24*3600*1000 = 2629746000 milliseconds)
timeSteps=c(2629746000)

r <- getRasterSlice(requestedParameter = "elevation",
                    lon_min = -13,
                    lon_max = 10,
                    lat_min = 40,
                    lat_max = 60,
                    requestedTimeSteps = NA,
                    date = "2020-01-01",
                    stacCatalogue = EDITOSTAC)

# with predefined user input
r <- getRasterSlice(requestedParameter = "elevation",
                    lon_min = -13,
                    lon_max = 10,
                    lat_min = 40,
                    lat_max = 60,
                    requestedTimeSteps = NA,
                    date = "2020-01-01",
                    stacCatalogue = EDITOSTAC,
                    select_layers = 1)

# with predefined user input
r <- getRasterSlice(requestedParameter = "Substrate",
                    lon_min = -13,
                    lon_max = 10,
                    lat_min = 40,
                    lat_max = 60,
                    requestedTimeSteps = NA,
                    stacCatalogue = EDITOSTAC)

plot(r)


