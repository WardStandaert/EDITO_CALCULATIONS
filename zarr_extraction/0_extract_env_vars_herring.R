setwd("~/work/EDITO_CALCULATIONS/zarr_extraction/")

# load functions, but also cached stac catalog named stacCatalog
source("editoTools.R")
options("outputdebug"=c('L','M'))

#the cached stacCatalog is called 'EDITOSTAC'
load(file = "./data-raw/editostacv2.par")

# the file to process
# datafile = "./data-raw/herring_test_set.csv"
datafile = "./data-raw/extract_test.csv"
pts = read.delim(datafile, sep=",") %>% sample_n(100)
glimpse(pts)

pts <- pts %>%
  mutate(Latitude = lat,
         Longitude = lon,
         Time = as.POSIXct(paste(1,pts$Month,pts$Year, sep = "-"), format = "%d-%m-%Y"))

glimpse(pts)

#the requested timestep resolution of the dataset in milliseconds
#in this case we work with monthly data (1 month = 30.436875*24*3600*1000 = 2629746000 milliseconds)
timeSteps=c(2629746000)


#the requested parameters, names in the stac catalogue table field 'par' .. see unique(stacCatalogue$par) for a list 

nms <- names(pts %>% select(-"...1",-X, -pa, -lon, -lat, -Longitude, -Latitude, -Year, -Month, -Time))
nms

parameters <- c("thetao", "so", "zooc", "phyc", "Energy", "Substrate", "elevation")


#check if they all exist
for ( parameter in parameters) {
  if(! parameter %in% unique(EDITOSTAC$par) ) 
    { dbl("Unknown parameter ", parameter)
  }
}

#add verbose= anything to get additional info on the positions ( par_x, par_y, par_z ) and time (par_t) found in the zarr files
enhanced_DF=enhanceDF(inputPoints = pts, 
                      requestedParameters = parameters, 
                      requestedTimeSteps = timeSteps, 
                      stacCatalogue = EDITOSTAC, 
                      verbose="on")

glimpse(enhanced_DF %>% select(all_of(c(nms,parameters)[order(c(nms,parameters))])))

par(mfrow = c(3,2))
plot(enhanced_DF$SST, enhanced_DF$thetao, cex.lab=1.3,
     xlim = c(min(enhanced_DF$SST, enhanced_DF$thetao), max(enhanced_DF$SST, enhanced_DF$thetao)),
     ylim = c(min(enhanced_DF$SST, enhanced_DF$thetao), max(enhanced_DF$SST, enhanced_DF$thetao)))
plot(enhanced_DF$SSS, enhanced_DF$so, cex.lab=1.3,
     xlim = c(min(enhanced_DF$SSS, enhanced_DF$so), max(enhanced_DF$SSS, enhanced_DF$so)),
     ylim = c(min(enhanced_DF$SSS, enhanced_DF$so), max(enhanced_DF$SSS, enhanced_DF$so)))
plot(enhanced_DF$Phyto, enhanced_DF$phyc, cex.lab=1.3,
     xlim = c(min(enhanced_DF$Phyto, enhanced_DF$phyc), max(enhanced_DF$Phyto, enhanced_DF$phyc)),
     ylim = c(min(enhanced_DF$Phyto, enhanced_DF$phyc), max(enhanced_DF$Phyto, enhanced_DF$phyc)))
plot(enhanced_DF$ZooPl, enhanced_DF$zooc, cex.lab=1.3,
     xlim = c(min(enhanced_DF$ZooPl, enhanced_DF$zooc), max(enhanced_DF$ZooPl, enhanced_DF$zooc)),
     ylim = c(min(enhanced_DF$ZooPl, enhanced_DF$zooc), max(enhanced_DF$ZooPl, enhanced_DF$zooc)))
plot(enhanced_DF$depth, enhanced_DF$elevation, cex.lab=1.3,
     xlim = c(min(enhanced_DF$depth, enhanced_DF$elevation), max(enhanced_DF$depth, enhanced_DF$elevation)),
     ylim = c(min(enhanced_DF$depth, enhanced_DF$elevation), max(enhanced_DF$depth, enhanced_DF$elevation)))

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

# write.csv(enhanced_DF, "tst/extract_test.csv")
glimpse(enhanced_DF)

sum(enhanced_DF$Energy_Description == enhanced_DF$ene_char) / length(enhanced_DF$Energy_Description)
sum(enhanced_DF$Substrate_Description == enhanced_DF$sub_char) / length(enhanced_DF$Substrate_Description)

par(mfrow = c(1,2))
plot(enhanced_DF$Energy, enhanced_DF$seabed_energy)
plot(enhanced_DF$Substrate, enhanced_DF$seabed_substrate)

enhanced_DF %>%
  group_by(Substrate_Description) %>% 
  count() %>% 
  full_join(enhanced_DF %>%
              group_by(sub_char) %>% 
              count(), by = c("Substrate_Description" = "sub_char")) %>%
  mutate(old_extraction = n.y,
         EDITO_extraction = n.x) %>%
  select(-n.x, -n.y)

enhanced_DF %>%
  group_by(Energy_Description) %>% 
  count() %>% 
  full_join(enhanced_DF %>%
              group_by(ene_char) %>% 
              count(), by = c("Energy_Description" = "ene_char")) %>%
  mutate(old_extraction = n.y,
         EDITO_extraction = n.x) %>%
  select(-n.x, -n.y)
