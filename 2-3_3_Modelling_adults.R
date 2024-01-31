library(tidyr)
library(dismo)
library(dplyr)
library(sf)
library(stringr)
library(mapview)
library(ape)
library(ggplot2)
library(tictoc)
library(RColorBrewer)
library(ENMeval)
library(tibble)
library(maxnet)
library(spdep)
library(dynamicSDM)
library(rJava)
library(rasterVis)
library(StatMatch)


#0. load data & functions ----

load("~/INPUT/final_st_lists.Rdata")
load("~/INPUT/bio_dat_co.Rdata")
load("~/INPUT/presences_after_spatial_filtering.rData")

setwd("~/OUTPUT/3.sp&env_filtered_model/adult/final_model_400")

all_vars <- str_remove_all(names(st_list_NEA2[[1]][[1]]) , "_\\d{4}_\\d{2}|_\\d{4}_\\d")
v_list_all <- all_vars[-which(all_vars == "Chl")]

#remove Phyto from stack list
st_list_NEA_cl <- st_list_NEA2
for (y in 1:length(2000:2020)) {
  for (m in 1:12) {
    st_list_NEA_cl[[y]][[m]] <- st_list_NEA2[[y]][[m]][[which(str_detect(names(st_list_NEA2[[y]][[m]]), paste(v_list_all, collapse = "|")))]]
  }
}

# Predictor variables
species <- tibble(scientific = c("Clupea harengus", "Scomber scombrus" , "Alosa fallax", "Dicentrarchus labrax"),
                  simple = c("herring", "mackerel" , "twaite_shad", "seabass"))

# Variables selected per species
# Twaite shad includes all vars except for phytoplankton because the variable is correlated with seabed_energy

v_list <- list()
v_list$herring <- c("SST", "SSS", "windfarms", "ZooPl", "Phyto", "EuphD", "seabed_energy", "seabed_substrate", "depth")
v_list$mackerel <- c("SST", "SSS", "windfarms", "max_SSV", "EuphD", "ZooPl", "depth")
v_list$twaite_shad <- c("Quarter_sin", "Quarter_cos", "SST", "SSS", "windfarms", "ZooPl", "seabed_energy", "seabed_substrate", "depth") 
v_list$seabass <- c("SST", "SSS", "windfarms", "ZooPl", "depth")


# Circular encoding for month where january (1) and december (12) are similar
reduced_dat <- reduced_dat %>%
  mutate(Month_sin = sin(Month * (2*pi/12)),
         Month_cos = cos(Month * (2*pi/12)),
         quarter = ceiling(Month / 3),
         Quarter_sin = sin(quarter * (2*pi/4)),
         Quarter_cos = cos(quarter * (2*pi/4)))


month_lvl <- tibble(mon_char = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                    Month = c(1:12))
substr_lvl <- tibble(sub_char = c("Fine mud", "Sand", "Muddy sand", "Mixed sediment",
                                  "Coarse substrate","Sandy mud or Muddy sand", "Seabed",
                                  "Rock or other hard substrata","Sandy mud", "Sandy mud or Muddy sand ",
                                  "Sediment","Fine mud or Sandy mud or Muddy sand"),
                     seabed_substrate = c(1:12))
energy_lvl <- tibble(ene_char = c("High energy", "Moderate energy", "Low energy", "No energy information"),
                     seabed_energy = c(1:4))

features_fun <- function(string) {
  all_options <- c("L","Q","H","P","T")
  key <- c("linear","quadratic","hinge","product","threshold")
  requested_options <- str_split(string, "")[[1]]
  ind1 <- which(all_options %in% requested_options)
  ind2 <- which(!all_options %in% requested_options)
  out <- c(paste0(key[ind1], "=true"), paste0(key[ind2], "=false"))
  return(out)
}

#0.1. environmental filtering ----
num_pr_vector <- c(400,400,160,160)

red_dat <- list()
for (s in 1:nrow(species)) {
  name <- species$simple[s]
  
  tmp_full_dat <- left_join(presences_list[[which(names(presences_list) == species$simple[s])]] %>% 
                              mutate(lon = Longitude, lat = Latitude), full_dat, by = c("lon","lat")) %>%
    filter(scientificname == species$scientific[s])
  
  df_pr <- tmp_full_dat[0,]
  df_rem <- tmp_full_dat
  
  #remove variables that have only one unique value
  tmp_df <- tmp_full_dat %>% 
    select(any_of(v_list[[s]])) %>%
    select(where(~ length(unique(.)) != 1))
  
  mal_dis <- mahalanobis.dist(tmp_df)
  
  # Select the two points that are furthest apart based on this distance and add to df_pr
  # (if plural pairs are present one pair is selected randomly)
  id <- which(mal_dis == max(mal_dis), arr.ind = TRUE)
  pnt_coord <- id[sample(nrow(id),1,replace=FALSE),]
  df_pr <- rbind(df_pr,df_rem[as.vector(pnt_coord)[1],], df_rem[as.vector(pnt_coord)[2],])
  # Remove this pair from the remaining points
  df_rem <- df_rem[-c(as.vector(pnt_coord)[1],as.vector(pnt_coord)[2]),]
  
  # This loop now selects the point that is furthest away from the current points in df_pr
  # in order to get a total of 75 points
  while (nrow(df_pr) < num_pr_vector[s]) {
    mal_dis <- mahalanobis.dist(df_pr[,12:16],df_rem[,12:16])
    mal_dis_tot <- colSums(mal_dis)
    id <- as.vector(which(mal_dis_tot == max(mal_dis_tot), useNames = FALSE))
    if (length(id) > 1) {id <- sample(id,1)}
    df_pr <- rbind(df_pr,df_rem[id,])
    df_rem <- df_rem[-id,]
  }
  red_dat[[s]] <- df_pr
  print(name)
}

reduced_dat <- bind_rows(red_dat)
table(reduced_dat$scientificname)
reduced_dat <- reduced_dat %>%
  mutate(Month_sin = sin(Month * (2*pi/12)),
         Month_cos = cos(Month * (2*pi/12)),
         quarter = ceiling(Month / 3),
         Quarter_sin = sin(quarter * (2*pi/4)),
         Quarter_cos = cos(quarter * (2*pi/4)))

#0.2. no environmental filtering ----
red_dat <- list()
for (s in 1:nrow(species)) {
  name <- species$simple[s]
  
  tmp_full_dat <- left_join(presences_list[[which(names(presences_list) == species$simple[s])]] %>% 
                              mutate(lon = Longitude, lat = Latitude), full_dat, by = c("lon","lat")) %>%
    filter(scientificname == species$scientific[s])
  
  red_dat[[s]] <- tmp_full_dat
  print(name)
}

reduced_dat <- bind_rows(red_dat)
table(reduced_dat$scientificname)

#1. extract environmental values for presences ---- 

# reduced_dat <- reduced_dat %>%
#   mutate(Month = as.factor(Month)) %>%
#   mutate(seabed_energy = as.factor(seabed_energy)) %>%
#   mutate(seabed_substrate = as.factor(seabed_substrate))

# full_dat <- full_dat %>%
#   mutate(Month = as.factor(Month)) %>%
#   mutate(seabed_energy = as.factor(seabed_energy)) %>%
#   mutate(seabed_substrate = as.factor(seabed_substrate))

#convert month to factor

pr_pa <- list()
pr_coord <- list()
pr_predv <- list()
pr_pa_coord_predv <- list()

for (s in 1:nrow(species)) {
  tmp_pres  <- reduced_dat %>%
    filter(scientificname == species$scientific[s])
  
  pr_pa[[s]] <- rep(1, nrow(tmp_pres))
  pr_coord[[s]] <- tmp_pres %>% dplyr::select(lon, lat, Year, Month)
  pr_predv[[s]] <- tmp_pres %>% dplyr::select(all_of(v_list[[which(names(v_list) == species$simple[s])]]))
  
  names(pr_pa)[s] <- names(pr_coord)[s] <- names(pr_predv)[s] <- species$simple[s]
}

lapply(pr_pa, head, n=2)              #pr points presence/absence values per species
lapply(pr_coord, head, n=2)           #pr points coordinates per species
lapply(pr_predv, head, n=2)           #pr points covariate extracts per species


#2.1. create background points and extract environmental values ----
set.seed(101)

bg_pa <- list()
bg_coord <- list()
bg_predv <- list()

for (s in 1:nrow(species)) {
  #create background points within ices regions that species occurs in
  tmp_dat <- reduced_dat %>%
    filter(scientificname == species$scientific[s])
  
  tmp_areas <- pull(unique(tmp_dat %>% dplyr::select(Area_27)))
  
  tmp_crop <- ices_shp[which(ices_shp$Area_27 %in% tmp_areas),]
  st_crs(tmp_crop) <- 3857
  tmp_crop <- st_transform(tmp_crop, st_crs("epsg:4326"))
  tmp_mask <- mask(final_common_mask, tmp_crop)
  crs(tmp_mask) <- 4326
  
  tmp_coords <- as.data.frame(coordinates(tmp_mask)[which(values(!is.na(tmp_mask))),])
  colnames(tmp_coords) <- c("lon","lat")
  
  tmp_background <- data.frame(tmp_coords[sample(1:nrow(tmp_coords), 10*length(pr_pa[[s]]), replace = TRUE),],
                               Year = sample(tmp_dat$Year, 10*length(pr_pa[[s]]), replace = TRUE),
                               Month = sample(unique(tmp_dat$Month), 10*length(pr_pa[[s]]), replace = TRUE)) %>%
    mutate(Month_sin = sin(Month * (2*pi/12)),
           Month_cos = cos(Month * (2*pi/12)),
           quarter = ceiling(Month / 3),
           Quarter_sin = sin(quarter * (2*pi/4)),
           Quarter_cos = cos(quarter * (2*pi/4)))

  #sample environmental layers
  tmp_bckgr <- data.frame()
  for (y in 1:length(2000:2020)) {
    for (m in 1:12) {
      tmp_coords_dat <- tmp_background %>%
        filter(Year == c(2000:2020)[y],
               Month == m)
      tmp_r <- st_list_NEA_cl[[y]][[m]]
      tmp_extract <- raster::extract(st_list_NEA_cl[[y]][[m]], 
                                     as.data.frame(tmp_coords_dat[,c(1,2)]), method = "simple")
      colnames(tmp_extract) <- colnames(tmp_extract) %>%
        str_remove_all("_\\d{4}_\\d{2}|_\\d{4}_\\d")
      tmp_df <- data.frame(tmp_coords_dat, tmp_extract)
      tmp_bckgr <- rbind(tmp_bckgr, tmp_df)
    }
    print(y)
  }
  
  #outputs
  bg_pa[[s]] <- rep(0, nrow(tmp_bckgr))
  
  bg_coord[[s]] <- tmp_bckgr %>% dplyr::select(lon, lat, Year, Month)
  
  bg_predv[[s]] <- tmp_bckgr %>%  #make seabed habitat a factor
    dplyr::select(all_of(v_list[[which(names(v_list) == species$simple[s])]])) #variable selection
  
   names(bg_pa)[s] <- names(bg_coord)[s] <- names(bg_predv)[s] <- species$simple[s]
}

lapply(bg_pa, head, n=2)              #bg points presence/absence values per species
lapply(bg_coord, head, n=2)           #bg points coordinates per species
lapply(bg_predv, head, n=2)           #bg points covariate extracts per species

#remove all tmp files
rm(list = ls()[which(str_detect(ls(), pattern = "^tmp_"))])


#visualize generated background points versus occurrences
mapview(bg_coord$herring[,1], bg_coord$herring[,2], crs = 'epsg:4326', col.regions = "red", layer.name = "background", cex = 3) +
  mapview(pr_coord$herring[,1], pr_coord$herring[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)
mapview(bg_coord$mackerel[,1], bg_coord$mackerel[,2], crs = 'epsg:4326', col.regions = "red", layer.name = "background", cex = 3) +
  mapview(pr_coord$mackerel[,1], pr_coord$mackerel[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)
mapview(bg_coord$twaite_shad[,1], bg_coord$twaite_shad[,2], crs = 'epsg:4326', col.regions = "red", layer.name = "background", cex = 3) +
  mapview(pr_coord$twaite_shad[,1], pr_coord$twaite_shad[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)
mapview(bg_coord$seabass[,1], bg_coord$seabass[,2], crs = 'epsg:4326', col.regions = "red", layer.name = "background", cex = 3) +
  mapview(pr_coord$seabass[,1], pr_coord$seabass[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)
  

mapview(pr_coord$herring[,1], pr_coord$herring[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)
mapview(pr_coord$mackerel[,1], pr_coord$mackerel[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)
mapview(pr_coord$twaite_shad[,1], pr_coord$twaite_shad[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)
mapview(pr_coord$seabass[,1], pr_coord$seabass[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)

# setwd("/data/home/innovauth/ward.standaert@vliz.be/OUTPUT/pa_coordinates")
# for (s in 1:nrow(species)) write.csv(pr_coord[[s]], paste0("adult_presences_",species$simple[s],".csv"))
# for (s in 1:nrow(species)) write.csv(bg_coord[[s]], paste0("adult_background_points_",species$simple[s],".csv"))

#2.2. create background points using maxdist of presences and extract environmental values ----
# set.seed(101)
# 
# bg_pa <- list()
# bg_coord <- list()
# bg_predv <- list()
# bg_pa_coord_predv <- list()
# 
# for (s in 1:nrow(species)) {
#   #create background points within ices regions that species occurs in
#   tmp_dat <- reduced_dat %>%
#     filter(scientificname == species$scientific[s])
#   
#   pts <- tmp_dat %>% dplyr::select(lon, lat)
#   coordinates(pts) <- ~lon+lat
#   projection(pts) <- CRS("+proj=longlat +datum=WGS84")
#   pts <- spTransform(pts, crs("+init=epsg:3034"))
#   
#   x <- polygons(circles(pts, d = 50000, lonlat = FALSE))
#   crs(x) <- crs("+init=epsg:3034")
#   x <- spTransform(x, crs("+proj=longlat +datum=WGS84"))
#   crs(final_common_mask) <- crs("+proj=longlat +datum=WGS84")
#   
#   tmp_mask <- mask(final_common_mask, x)
#   tmp_pts <- randomPoints(tmp_mask, 10*length(pr_pa[[s]]))
#   
#   tmp_background <- data.frame(lon = tmp_pts[,1],
#                                lat = tmp_pts[,2],
#                                Year = sample(c(2000:2020), 10*length(pr_pa[[s]]), replace = T),
#                                Month = sample(unique(tmp_dat$Month), 10*length(pr_pa[[s]]), replace = T))
#   
#   #sample environmental layers
#   tmp_bckgr <- data.frame()
#   for (y in 1:length(2000:2020)) {
#     for (m in 1:12) {
#       tmp_coords_dat <- tmp_background %>%
#         filter(Year == c(2000:2020)[y],
#                Month == m)
#       tmp_r <- st_list_NEA_cl[[y]][[m]]
#       tmp_extract <- raster::extract(st_list_NEA_cl[[y]][[m]], 
#                                      as.data.frame(tmp_coords_dat[,c(1,2)]), method = "simple")
#       colnames(tmp_extract) <- colnames(tmp_extract) %>%
#         str_remove_all("_\\d{4}_\\d{2}|_\\d{4}_\\d")
#       tmp_df <- data.frame(tmp_coords_dat, tmp_extract)
#       tmp_bckgr <- rbind(tmp_bckgr, tmp_df)
#     }
#     print(y)
#   }
#   
#   #outputs
#   bg_pa[[s]] <- rep(0, nrow(tmp_bckgr))
#   
#   bg_coord[[s]] <- tmp_bckgr %>% dplyr::select(lon, lat, Year, Month)
#   
#   bg_predv[[s]] <- tmp_bckgr %>% 
#     mutate(Month = as.factor(Month),      #make month a factor
#            seabed_energy = as.factor(seabed_energy),      #make seabed energy a factor
#            seabed_substrate = as.factor(seabed_substrate)) %>%  #make seabed habitat a factor
#     dplyr::select(all_of(v_list[[which(names(v_list) == species$simple[s])]])) #variable selection
#   
#   bg_pa_coord_predv[[s]] <- tmp_bckgr %>% 
#     mutate(pa = rep(0, nrow(tmp_bckgr)),
#            Month = as.factor(Month)) %>%  #make month a factor
#     dplyr::select(all_of(c("pa", "lon", "lat", v_list[[which(names(v_list) == species$simple[s])]]))) 
#   
#   names(bg_pa)[s] <- names(bg_coord)[s] <- names(bg_predv)[s] <- names(bg_pa_coord_predv)[s] <- species$simple[s]
# }
# 
# lapply(bg_pa, head, n=2)              #bg points presence/absence values per species
# lapply(bg_coord, head, n=2)           #bg points coordinates per species
# lapply(bg_predv, head, n=2)           #bg points covariate extracts per species
# lapply(bg_pa_coord_predv, head, n=2)  #bg points full combination per species
# 
# rm(tmp_extract, tmp_coords_dat, tmp_r, tmp_df, tmp_bckgr, 
#    tmp_areas, tmp_crop, tmp_pts, tmp_mask, tmp_dat, tmp_background)
# 
# 
# #visualize generated background points versus occurrences
# mapview(bg_coord$herring[,1], bg_coord$herring[,2], crs = 'epsg:4326', col.regions = "red", layer.name = "background", cex = 3) +
#   mapview(pr_coord$herring[,1], pr_coord$herring[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)
# mapview(bg_coord$mackerel[,1], bg_coord$mackerel[,2], crs = 'epsg:4326', col.regions = "red", layer.name = "background", cex = 3) +
#   mapview(pr_coord$mackerel[,1], pr_coord$mackerel[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)
# mapview(bg_coord$twaite_shad[,1], bg_coord$twaite_shad[,2], crs = 'epsg:4326', col.regions = "red", layer.name = "background", cex = 3) +
#   mapview(pr_coord$twaite_shad[,1], pr_coord$twaite_shad[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)
# mapview(bg_coord$seabass[,1], bg_coord$seabass[,2], crs = 'epsg:4326', col.regions = "red", layer.name = "background", cex = 3) +
#   mapview(pr_coord$seabass[,1], pr_coord$seabass[,2], crs = 'epsg:4326', col.regions = "green", layer.name = "occurrences", cex = 3)
# 
# r <- st_list_NEA_cl[[1]][[1]][[1]]
# crs(r) <- crs("+proj=longlat +datum=WGS84")
# mapview(r) +
#   mapview(bg_coord$herring[,1], bg_coord$herring[,2], crs = 'epsg:4326', col.regions = "red", layer.name = "background", cex = 3)
#   

# start from SAVE
save(pr_pa, pr_coord, pr_predv,
     bg_pa, bg_coord, bg_predv,
     file = "SAVE/NEA_presence_background.Rdata")

# SAVE ----
load("SAVE/NEA_presence_background.Rdata")

#3 ENMeval model creation  ----
#convert month to character bc algorithm cant handle factors that look like numeric
eval_res_list <- list()
tmp_pr_predv <- list()
tmp_bg_predv <- list()
for (s in 1:nrow(species)) {
  if(all(c("seabed_energy","seabed_substrate","Month") %in% names(pr_predv[[s]]))) {
    tmp_pr_predv[[s]] <- pr_predv[[s]] %>%
      left_join(month_lvl, by = "Month") %>%
      left_join(energy_lvl, by = "seabed_energy") %>%
      left_join(substr_lvl, by = "seabed_substrate") %>%
      dplyr::select(-Month, -seabed_energy, -seabed_substrate)
    
    tmp_bg_predv[[s]] <- bg_predv[[s]] %>%
      left_join(month_lvl, by = "Month") %>%
      left_join(energy_lvl, by = "seabed_energy") %>%
      left_join(substr_lvl, by = "seabed_substrate") %>%
      dplyr::select(-Month, -seabed_energy, -seabed_substrate)
    cat <- c("mon_char", "ene_char","sub_char")
      }
  else if(all(c("seabed_energy","seabed_substrate") %in% names(pr_predv[[s]]))) {
    tmp_pr_predv[[s]] <- pr_predv[[s]] %>%
      left_join(energy_lvl, by = "seabed_energy") %>%
      left_join(substr_lvl, by = "seabed_substrate") %>%
      dplyr::select(-seabed_energy, -seabed_substrate)
    
    tmp_bg_predv[[s]] <- bg_predv[[s]] %>%
      left_join(energy_lvl, by = "seabed_energy") %>%
      left_join(substr_lvl, by = "seabed_substrate") %>%
      dplyr::select(-seabed_energy, -seabed_substrate)
    cat <- c("ene_char","sub_char")
  }
  else if("Month" %in% names(pr_predv[[s]])) {
    tmp_pr_predv[[s]] <- pr_predv[[s]] %>%
      left_join(month_lvl, by = "Month") %>%
      dplyr::select(-Month)
    
    tmp_bg_predv[[s]] <- bg_predv[[s]] %>%
      left_join(month_lvl, by = "Month") %>%
      dplyr::select(-Month)
    cat <- "mon_char"
  }
  else {
    tmp_pr_predv[[s]] <- pr_predv[[s]]
    tmp_bg_predv[[s]] <- bg_predv[[s]]
    cat <- NULL
    
  }

  #if any column exists with only one unique value --> delete it
  if(any(apply(rbind(tmp_pr_predv[[s]], tmp_bg_predv[[s]]), 2, function(x) length(unique(x)) == 1))) {
    tmp_ind <- which(apply(rbind(tmp_pr_predv[[s]], tmp_bg_predv[[s]]), 2, function(x) length(unique(x)) == 1))
    tmp_ind_v <- names(tmp_ind)
    tmp_pr_predv[[s]] <- tmp_pr_predv[[s]][-tmp_ind]
    tmp_bg_predv[[s]] <- tmp_bg_predv[[s]][-tmp_ind]
    v_list[[s]] <- v_list[[s]][-which(v_list[[s]] == tmp_ind_v)]
    print(paste0("removed ", names(tmp_ind), " as a variable for ", species$simple[s]))
  }
  
  tmp_pr_predv[[s]] <- data.frame(pr_coord[[s]] %>% dplyr::select(lon, lat), tmp_pr_predv[[s]]) #input for ENMevaluate requires lon & lat as first covariates
  tmp_bg_predv[[s]] <- data.frame(bg_coord[[s]] %>% dplyr::select(lon, lat), tmp_bg_predv[[s]]) #input for ENMevaluate requires lon & lat as first covariates
  
  eval_res_list[[s]] <- ENMeval::ENMevaluate(occs = tmp_pr_predv[[s]],
                                    bg = tmp_bg_predv[[s]],
                                    tune.args = list(fc = c("L","LQ","LQH"),
                                                     rm = c(1,2,4,8,32)),
                                    algorithm = "maxnet",
                                    partitions = "randomkfold",
                                    categoricals = cat,
                                    doClamp = TRUE,
                                    parallel = TRUE)
  
  
  print(paste0(species$simple[s], " done"))
}

save(eval_res_list, file = "SAVE/NEA_ENMeval_outcomes.Rdata")

# SAVE ----
load("SAVE/NEA_ENMeval_outcomes.Rdata")

for(s in 1:nrow(species)) print(eval_res_list[[s]] %>%  eval.results() %>% filter(delta.AICc == 0))


# ## for full dataset -----
# tmp_pr_predv <- list()
# tmp_bg_predv <- list()
# model_list <- list()
# for (s in 1:nrow(species)) {
#   
#     tmp_pr_predv[[s]] <- pr_predv[[s]]
#     tmp_bg_predv[[s]] <- bg_predv[[s]]
#   
#   #if any column exists with only one unique value --> delete it
#   if(any(apply(rbind(tmp_pr_predv[[s]], tmp_bg_predv[[s]]), 2, function(x) length(unique(x)) == 1))) {
#     tmp_ind <- which(apply(rbind(tmp_pr_predv[[s]], tmp_bg_predv[[s]]), 2, function(x) length(unique(x)) == 1))
#     tmp_ind_v <- names(tmp_ind)
#     tmp_pr_predv[[s]] <- tmp_pr_predv[[s]][-tmp_ind]
#     tmp_bg_predv[[s]] <- tmp_bg_predv[[s]][-tmp_ind]
#     v_list[[s]] <- v_list[[s]][-which(v_list[[s]] == tmp_ind_v)]
#     print(paste0("removed ", names(tmp_ind), " as a variable for ", species$simple[s]))
#   }
#     
#   if(all(c("ene_char", "sub_char") %in% names(tmp_pr_predv[[s]]))) x <- rbind(tmp_pr_predv[[s]],tmp_bg_predv[[s]]) %>%
#       mutate(ene_char = as.numeric(ene_char),
#              sub_char = as.numeric(sub_char))
#       
#   else x <- rbind(tmp_pr_predv[[s]],tmp_bg_predv[[s]])
# 
#   eval_res <- eval_res_list[[s]]
#   res <- eval.results(eval_res)
#   opt.aicc <- res %>% filter(delta.AICc == 0)
# 
#   rm <- levels(opt.aicc$rm)[as.numeric(opt.aicc$rm[1])]
#   f <- features_fun(opt.aicc$fc)
#   
#   
#   model_list[[s]] <- dismo::maxent(x,
#                             p = c(pr_pa[[s]], bg_pa[[s]]),
#                             factors = c("ene_char","sub_char"),
#                             args = c(f,
#                                      paste0("beta_threshold=", rm), paste0("beta_categorical=", rm), 
#                                      paste0("beta_lqp=", rm), paste0("beta_hinge=", rm)),
#                             path = "/data/home/innovauth/ward.standaert@vliz.be/ENMevaluate outcomes/spatial and environmental filtering")
#     
#   names(model_list)[[s]] <- species$simple[s]
#   print(paste("finished:", species$simple[s]))
# }




##gives different results --> try old method to do jackknife (thesis)
# model_list_enm <- list()
# for (s in 1:nrow(species)) {
#   eval_res <- eval_res_list[[s]]
#   res <- eval.results(eval_res)
#   opt.aicc <- res %>% filter(delta.AICc == 0)
#   
#   model_list_enm[[s]]      <- maxent(x = rbind(pr_predv[[s]],bg_predv[[s]]),
#                                      p = c(pr_pa[[s]], bg_pa[[s]]),
#                                      args = c("jackknife=true", 
#                                               features_fun(opt.aicc$fc), #optimal hinges
#                                               paste0("beta_hinge=", as.numeric(opt.aicc$rm))))
#   names(model_list_enm)[[s]] <- species$simple[s]
#   print(paste("finished:", species$simple[s]))
# }
# 
# plot(model_list_enm[[1]])
# response(model_list_enm[[1]])
# 
# r1 <- st_list_NEA_cl[[1]][[1]][[1]]
# values(r1) <- 8
# names(r1) <- "mon_char"
# r2 <- st_list_NEA_cl[[1]][[1]][[1]]
# values(r2) <- 2000
# names(r2) <- "Year"
# plot_st <- stack(st_list_NEA_cl[[1]][[1]], r1, r2)
# names(plot_st) <- str_remove(names(plot_st), "_\\d{4}_\\d")
# names(plot_st)[which(names(plot_st) == "seabed_energy" | names(plot_st) == "seabed_substrate")] <- c("ene_char","sub_char")
# 
# # reclass_table_energy <- energy_lvl$seabed_energy
# # names(energy_lvl) <- c("to","from")
# # qualitative_raster <- reclassify(st_list_NEA_cl[[1]][[1]][[7]], cbind(energy_lvl$seabed_energy, energy_lvl$ene_char))
# 
# prediction1 <- predict(plot_st , eval_res@models$fc.L_rm.1, type = "cloglog", 
#                        factors = list(mon_char=factor(month_lvl$mon_char, levels = month_lvl$mon_char),
#                                       ene_char=factor(substr_lvl$sub_char, levels = substr_lvl$sub_char),
#                                       sub_char=factor(energy_lvl$ene_char, levels = energy_lvl$ene_char)))
# prediction1 <- predict(plot_st , eval_res@models$fc.L_rm.1, type = "cloglog")
# plot(prediction1)
# 
# 
# 
# factor(substr_lvl$sub_char, levels = substr_lvl$sub_char)


#3.2 SAC corrected ENMeval model creation  ----
# #doesnt work
# #convert month to character bc algorithm cant handle factors that look like numeric
# month_lvl <- data.frame(mon_char = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
#                         Month = factor(c(1:12)))
# substr_lvl <- data.frame(sub_char = c("Fine mud", "Sand", "Muddy sand",
#                                       "Mixed sediment","Coarse substrate","Sandy mud or Muddy sand",
#                                       "Seabed","Rock or other hard substrata","Sandy mud",
#                                       "Sandy mud or Muddy sand ","Sediment","Fine mud or Sandy mud or Muddy sand"),
#                          seabed_substrate = factor(c(1:12)))
# energy_lvl <- data.frame(ene_char = c("High energy", "Moderate energy", "Low energy", "No energy information"),
#                          seabed_energy = factor(c(1:4)))
# 
# 
# eval_res_list <- list()
# tmp_pr_predv <- list()
# tmp_bg_predv <- list()
# for (s in 1:nrow(species)) {
#   if(all(c("seabed_energy","seabed_substrate","Month") %in% names(pr_predv[[s]]))) {
#     tmp_pr_predv[[s]] <- pr_predv[[s]] %>%
#       left_join(month_lvl, by = "Month") %>%
#       left_join(energy_lvl, by = "seabed_energy") %>%
#       left_join(substr_lvl, by = "seabed_substrate") %>%
#       dplyr::select(-Month, -seabed_energy, -seabed_substrate)
#     
#     tmp_bg_predv[[s]] <- bg_predv[[s]] %>%
#       left_join(month_lvl, by = "Month") %>%
#       left_join(energy_lvl, by = "seabed_energy") %>%
#       left_join(substr_lvl, by = "seabed_substrate") %>%
#       dplyr::select(-Month, -seabed_energy, -seabed_substrate)
#   }
#   else if(all(c("seabed_energy","seabed_substrate") %in% names(pr_predv[[s]]))) {
#     tmp_pr_predv[[s]] <- pr_predv[[s]] %>%
#       left_join(energy_lvl, by = "seabed_energy") %>%
#       left_join(substr_lvl, by = "seabed_substrate") %>%
#       dplyr::select(-seabed_energy, -seabed_substrate)
#     
#     tmp_bg_predv[[s]] <- bg_predv[[s]] %>%
#       left_join(energy_lvl, by = "seabed_energy") %>%
#       left_join(substr_lvl, by = "seabed_substrate") %>%
#       dplyr::select(-seabed_energy, -seabed_substrate)
#   }
#   else if("Month" %in% names(pr_predv[[s]])) {
#     tmp_pr_predv[[s]] <- pr_predv[[s]] %>%
#       left_join(month_lvl, by = "Month") %>%
#       dplyr::select(-Month)
#     
#     tmp_bg_predv[[s]] <- bg_predv[[s]] %>%
#       left_join(month_lvl, by = "Month") %>%
#       dplyr::select(-Month)
#   }
#   else {
#     tmp_pr_predv[[s]] <- pr_predv[[s]]
#     tmp_bg_predv[[s]] <- bg_predv[[s]]
#   }
#   
#   #add autocovariate
#   ac <- autocov_dist(c(pr_pa[[s]], bg_pa[[s]]), 
#                      rbind(pr_coord[[s]][,c(1,2)], bg_coord[[s]][,c(1,2)]), 
#                      nbs = 1, type = "inverse", zero.policy = TRUE)
#   ac[which(ac == Inf)] <- 0
# 
#   tmp_pr_predv[[s]]$ac <- ac[1:100]
#   tmp_bg_predv[[s]]$ac <- ac[101:1100]
#   
#   eval_res_list[[s]] <- ENMevaluate(occs = tmp_pr_predv[[s]], 
#                                     bg = tmp_bg_predv[[s]],
#                                     tune.args = list(fc = c("L","LQ","LQH","H"), rm = c(1,2,4,8,32)),
#                                     algorithm = "maxnet",
#                                     partitions = "randomkfold",
#                                     doClamp = TRUE,
#                                     parallel = TRUE)
#   print(paste0(species$simple[s], " done"))
# }
# 
# setwd("/data/home/innovauth/ward.standaert@vliz.be/ENMevaluate outcomes/spatial and environmental filtering with maxdist/")
# save(eval_res_list, file = "ENMeval_spec_vars_d_no-m_L_LQ_LQH_H_1_5.Rdata")
# 
# setwd("/data/home/innovauth/ward.standaert@vliz.be/ENMevaluate outcomes/spatial and environmental filtering with maxdist/L,LQ,LQH,H - 1:5 spec vars w depth no month")
# write.csv(eval.results(eval_res_list[[1]]), "ENMeval_results_herring.csv")
# write.csv(eval.results(eval_res_list[[2]]), "ENMeval_results_mackerel.csv")
# write.csv(eval.results(eval_res_list[[3]]), "ENMeval_results_twaite_shad.csv")
# write.csv(eval.results(eval_res_list[[4]]), "ENMeval_results_seabass.csv")
# 
# for(s in 1:nrow(species)) print(eval_res_list[[s]] %>%  eval.results() %>% filter(delta.AICc == 0))


#4 model evaluation AUC & TSS ----

AUC_maxent <- list()
TSS_maxent <- list()

for (s in 1:nrow(species)) {
  tic(species$simple[s])
  Presences <- pr_predv[[s]]
  Background <- bg_predv[[s]]
  
  tmp_AUC <- list()
  tmp_TSS <- list()
  for (i in 1:10) {
    tic(paste("iteration", i))
    #We generate "k" groups:
    k<-4 # This value can be modified. If we choose a 4, our divison of the data will be 75% Vs 25%
    # there is no set rule about the k-groups. Arbitrary.
    
    groups_pres<-kfold(Presences,k) #Kfold divide the data, assigning every row to one of the K groups randomly. 
    groups_abs<-kfold(Background,k)
    
    #Four groups will be used to generate the model and the rest of the point (one group) will be used to evaluate it:
    EvalBg<-Background[groups_pres==1,]
    TrainBg<-Background[groups_pres!=1,]
    
    EvalPres<-Presences[groups_pres==1,]
    TrainPres<-Presences[groups_pres!=1,]

    #get model settings
    eval_res <- eval_res_list[[s]]
    res <- eval.results(eval_res)
    opt.aicc <- res %>% filter(delta.AICc == 0)

    rm <- levels(opt.aicc$rm)[as.numeric(opt.aicc$rm[1])]
    f <- features_fun(opt.aicc$fc)
    
    
    maxent_model <- maxent(x = rbind(TrainPres, TrainBg),
                                  p = c(rep(1, nrow(TrainPres)), rep(0, nrow(TrainBg))),
                                  factors = c("seabed_energy","seabed_substrate","Month"),
                                  args = c(f,
                                           paste0("beta_threshold=", rm), paste0("beta_categorical=", rm),
                                           paste0("beta_lqp=", rm), paste0("beta_hinge=", rm)))
    

    EvalBgRes  <- predict(maxent_model, EvalBg, type = "cloglog")
    EvalPresRes  <- predict(maxent_model, EvalPres, type = "cloglog")
    
    tmp_AUC[[i]]<-evaluate(c(EvalPresRes), c(EvalBgRes))
    umbral_maxent<-threshold( tmp_AUC[[i]])   #creemos que 'thershold' es una funci?n que calcula el umbral que maximiza el kappa
    tmp_TSS[[i]]<-evaluate(p=c(EvalPresRes), a=c(EvalBgRes), tr=umbral_maxent$spec_sens) # umbral_maxent$spec_sens: gives you the threshold to maximize the sum of the specificity and sensitivity
    
    toc()
  }
  
  AUC_maxent[[s]] <- tmp_AUC
  TSS_maxent[[s]] <- tmp_TSS
  names(AUC_maxent)[s] <- names(TSS_maxent)[s] <- species$simple[s]
  toc()
}

# #Lets see the results! 
AUC_maxent_vect   <- data.frame(matrix(nrow = 10, ncol = 4))
TPR_maxent_vect   <- data.frame(matrix(nrow = 10, ncol = 4))
TNR_maxent_vect   <- data.frame(matrix(nrow = 10, ncol = 4))

for (s in 1:nrow(species)) {
  AUC_maxent_vect[,s]  <- sapply(AUC_maxent[[s]],function(x){slot(x,'auc')})
  TPR_maxent_vect[,s]  <- sapply(TSS_maxent[[s]],function(x){slot(x,'TPR')})
  TNR_maxent_vect[,s]  <- sapply(TSS_maxent[[s]],function(x){slot(x,'TNR')})
  
  colnames(AUC_maxent_vect)[s] <- colnames(TPR_maxent_vect)[s] <- 
    colnames(TNR_maxent_vect)[s] <-  species$simple[s]
}
TSS_maxent <- TPR_maxent_vect + TNR_maxent_vect - 1

#AUC: 
# at 0.5 your model is no better than random (throwing the dice)
# at 1 your model is good
AUC_maxent_df <- AUC_maxent_vect %>%
  gather(key = "variable",
         value = "value")

plot(factor(AUC_maxent_df$variable, 
            levels = c("herring", "mackerel", "seabass", "twaite_shad"), 
            labels = c("Atlantic herring", "Atlantic mackerel", "European seabass", "Twaite shad")), 
     AUC_maxent_df$value, main = "AUC", xlab = "", ylab = "")

#True skill statistic
TSS_maxent_df <- TSS_maxent %>%
  gather(key = "variable",
         value = "value")

plot(factor(TSS_maxent_df$variable, 
            levels = c("herring", "mackerel", "seabass", "twaite_shad"), 
            labels = c("Atlantic herring", "Atlantic mackerel", "European seabass", "Twaite shad")), 
     TSS_maxent_df$value, main = "TSS", xlab = "", ylab = "")

sum_auc <- dplyr::group_by(AUC_maxent_df, variable) %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value))
sum_tss <- dplyr::group_by(TSS_maxent_df, variable) %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value))

TPR_maxent_df <- TPR_maxent_vect %>%
  gather(key = "variable",
         value = "value")

TNR_maxent_df <- TNR_maxent_vect %>%
  gather(key = "variable",
         value = "value")

sum_TPR <- dplyr::group_by(TPR_maxent_df, variable) %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value))
sum_TNR <- dplyr::group_by(TNR_maxent_df, variable) %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value))

#t.test(TPR_maxent_vect, TNR_maxent_vect)

sum_all <- cbind(sum_auc, sum_tss[,c(2:3)], sum_TPR[,c(2:3)], sum_TNR[,c(2:3)])
colnames(sum_all) <- c("species",
                       "av_AUC", "sd_AUC", 
                       "av_TSS", "sd_TSS",
                       "av_TPR", "sd_TPR", 
                       "av_TNR", "sd_TNR")
sum_all[match(sum_all$species, species$simple),]


write.csv(sum_all, "evaluation/NEA_validation_metrics.csv")

myColors <- c("#332288", "#DDCC77", "#CC6677", "grey33") ## this colorpalette is distinguishable also for people colorblindness (https://davidmathlogic.com/colorblind/#%23332288-%23117733-%2344AA99-%2388CCEE-%23DDCC77-%23CC6677-%23AA4499-%23882255)
plot(factor(AUC_maxent_df$variable, 
            levels = c("herring", "mackerel", "seabass", "twaite_shad"), 
            labels = c("Atlantic herring", "Atlantic mackerel", "European seabass", "Twaite shad")), 
     AUC_maxent_df$value, main = "AUC", xlab = "", ylab = "")

custom_labels <- c("Atlantic\nherring", "Atlantic\nmackerel", "European\nseabass", "Twaite\nshad")

## boxplots ####
png("evaluation/NEA_AUC_boxplot.png", width = 6, height =6, res = 100, units = 'in')
par(las = 2, mar = c(6, 5, 2, 2), mgp = c(3, 0.5, 0))
boxplot(
  AUC_maxent_df$value ~ factor(AUC_maxent_df$variable,
                               levels = c("herring", "mackerel", "seabass", "twaite_shad")),
  col = myColors,
  ylab = "AUC",
  xlab = "",
  cex.axis = 1.2,
  cex.lab = 1.2,
  names = custom_labels  # Use custom labels here
)
dev.off()

png("evaluation/NEA_TSS_boxplot.png", width = 6, height =6, res = 100, units = 'in')
par(las = 2, mar = c(6, 5, 2, 2), mgp = c(3, 0.5, 0))
boxplot(
  TSS_maxent_df$value ~ factor(TSS_maxent_df$variable,
                               levels = c("herring", "mackerel", "seabass", "twaite_shad")),
  col = myColors,
  ylab = "TSS",
  xlab = "",
  cex.axis = 1.2,
  cex.lab = 1.2,
  names = custom_labels  # Use custom labels here
)
dev.off()


#5. Spatial autocorrelation ----

## presence and background data ----
# -> how clumped is data compaired to background 
morans_I_list <- list()
for (s in 1:nrow(species)) {
  tic(paste0(species$simple[s], " - done"))
  
  eval_res <- eval_res_list[[s]]
  res <- eval.results(eval_res)
  opt.aicc <- res %>% filter(delta.AICc == 0)
  mod.seq <- eval.models(eval_res)[[opt.aicc$tune.args]]
  
  #2. make original prediction
  df_pred <- rbind(pr_predv[[s]], bg_predv[[s]])
  if(all(c("seabed_energy", "seabed_substrate", "Month") %in% names(df_pred))) {
    df_pred <- df_pred %>%
      mutate(sub_char = as.numeric(seabed_substrate)) %>%
      mutate(ene_char = as.numeric(seabed_energy)) %>%
      mutate(mon_char = as.numeric(Month)) %>%
      dplyr::select(-seabed_substrate, -seabed_energy, -Month)
  }
  if(all(c("seabed_energy", "seabed_substrate") %in% names(df_pred))) {
    df_pred <- df_pred %>%
      mutate(sub_char = as.numeric(seabed_substrate)) %>%
      mutate(ene_char = as.numeric(seabed_energy)) %>%
      dplyr::select(-seabed_substrate, -seabed_energy)
  }

  pred_vals <- predict(mod.seq, df_pred, se.fit=TRUE, type = "cloglog")
  
  res <- c(pr_pa[[s]],bg_pa[[s]]) - pred_vals
  
  dists <- as.matrix(dist(rbind(cbind(pr_coord[[s]][,1], pr_coord[[s]][,2]),
                                cbind(bg_coord[[s]][,1], bg_coord[[s]][,2]))))
  dists.inv <- 1/dists
  diag(dists.inv) <- 0
  dists.inv[is.infinite(dists.inv)] <- 0 #remove infinite values
  
  #Global Moran's I from ape package (use prob)
  #model residuals of maxent are the occurrence-probabilities predicted values (Mateo-Tomás & Olea 2010)
  morans_I_list[[s]] <- Moran.I(c(res), dists.inv, scaled = TRUE, alternative = "greater")
  toc()
}

#for full dataset ----
# morans_I_list <- list()
# for (s in 1:nrow(species)) {
#   tic(paste0(species$simple[s], " - done"))
# 
#   #2. make original prediction
#   df_pred <- rbind(pr_predv[[s]], bg_predv[[s]])
#   
#   pred_vals <- predict(model_list[[s]], df_pred, se.fit=TRUE, type = "cloglog")
#   
#   res <- c(pr_pa[[s]],bg_pa[[s]]) - pred_vals
#   
#   dists <- as.matrix(dist(rbind(cbind(pr_coord[[s]][,1], pr_coord[[s]][,2]),
#                                 cbind(bg_coord[[s]][,1], bg_coord[[s]][,2]))))
#   dists.inv <- 1/dists
#   diag(dists.inv) <- 0
#   dists.inv[is.infinite(dists.inv)] <- 0 #remove infinite values
#   
#   #Global Moran's I from ape package (use prob)
#   #model residuals of maxent are the occurrence-probabilities predicted values (Mateo-Tomás & Olea 2010)
#   morans_I_list[[s]] <- Moran.I(c(res), dists.inv, scaled = TRUE, alternative = "greater")
#   toc()
# }


morans_I_list[[1]]   
#After sp & geo filtering:  Moran's I value of 0.06, expected -0.0009; significant
morans_I_list[[2]]   
#After sp & geo filtering:  Moran's I value of 0.06, expected -0.0009; significant
morans_I_list[[3]]
#After sp & geo filtering:  Moran's I value of 0.03, expected -0.0009; significant
morans_I_list[[4]]
#After sp & geo filtering:  Moran's I value of 0.09, expected -0.0009; significant

x <- cbind(rbind(pr_coord[[s]], bg_coord[[s]]), pred_vals)
coordinates(x) <- ~lon+lat
mapview(x, zcol = "pred_vals")


## presence data ----
morans_I_list <- list()
for (s in 1:nrow(species)) {
  tic(paste0(species$simple[s], " - done"))
  eval_res <- eval_res_list[[s]]
  res <- eval.results(eval_res)
  opt.aicc <- res %>% filter(delta.AICc == 0)
  mod.seq <- eval.models(eval_res)[[opt.aicc$tune.args]]
  
  #2. make original prediction
  df_pred <- pr_predv[[s]]
  if(all(c("seabed_energy", "seabed_substrate", "Month") %in% names(df_pred))) {
    df_pred <- df_pred %>%
      left_join(substr_lvl, by = "seabed_substrate") %>%
      left_join(energy_lvl, by = "seabed_energy") %>%
      left_join(month_lvl, by = "Month") %>%
      dplyr::select(-seabed_substrate, -seabed_energy, -Month)
  }
  if(all(c("seabed_energy", "seabed_substrate") %in% names(df_pred))) {
      df_pred <- df_pred %>%
        mutate(sub_char = as.numeric(seabed_substrate)) %>%
        mutate(ene_char = as.numeric(seabed_energy)) %>%
        dplyr::select(-seabed_substrate, -seabed_energy)
    }
  if (c("Month") %in% names(df_pred)) df_pred <- df_pred %>%
    left_join(month_lvl, by = "Month") %>%
    dplyr::select(-Month)
  
  #get model
  eval_res <- eval_res_list[[s]]
  opt.aicc <- eval.results(eval_res_list[[s]]) %>% filter(delta.AICc == 0)
  mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]    
  
  pred_vals <- predict(mod, df_pred, se.fit=TRUE, type = "cloglog")
  
  res <- pr_pa[[s]] - pred_vals
  
  dists <- as.matrix(dist(cbind(pr_coord[[s]][,1], pr_coord[[s]][,2])))
  dists.inv <- 1/dists
  diag(dists.inv) <- 0
  dists.inv[is.infinite(dists.inv)] <- 0 #remove infinite values
  
  #Global Moran's I from ape package
  morans_I_list[[s]] <- Moran.I(c(pred_vals), dists.inv, scaled = TRUE, alternative = "greater")
  
  #Global Moran's I from spdep package - takes longer and gives same result
  # morans_I_list_spdep[[s]] <- moran.test(res,mat2listw(dists.inv, style = "W"))
  toc()
}


# #for full dataset ----
# morans_I_list <- list()
# for (s in 1:nrow(species)) {
#   tic(paste0(species$simple[s], " - done"))
#   
#   #2. make original prediction
#   df_pred <- pr_predv[[s]]
#   
#   pred_vals <- predict(model_list[[s]], df_pred, se.fit=TRUE, type = "cloglog")
#   
#   res <- pr_pa[[s]] - pred_vals
#   
#   dists <- as.matrix(dist(rbind(cbind(pr_coord[[s]][,1], pr_coord[[s]][,2]))))
#   dists.inv <- 1/dists
#   diag(dists.inv) <- 0
#   dists.inv[is.infinite(dists.inv)] <- 0 #remove infinite values
#   
#   #Global Moran's I from ape package (use prob)
#   #model residuals of maxent are the occurrence-probabilities predicted values (Mateo-Tomás & Olea 2010)
#   morans_I_list[[s]] <- Moran.I(c(res), dists.inv, scaled = TRUE, alternative = "greater")
#   toc()
# }


morans_I_list[[1]]   
#Before filtering:          Moran's I value of 0.18, expected -0.00008; significant
#After sp filtering:        Moran's I value of 0.11, expected -0.00008; significant
#After sp & geo filtering:  Moran's I value of 0.09, expected -0.01; significant
morans_I_list[[2]]   
#Before filtering:          Moran's I value of 0.16, expected -0.0001; significant
#After sp & geo filtering:  Moran's I value of 0.10, expected -0.01; significant
morans_I_list[[3]]
#Before filtering:          Moran's I value of 0.20, expected -0.002; significant
#After sp & geo filtering:  Moran's I value of 0.19, expected -0.01; significant
morans_I_list[[4]]
#Before filtering:          Moran's I value of 0.21, expected -0.001; significant
#After sp & geo filtering:  Moran's I value of 0.21, expected -0.01; significant

# setwd("/data/home/innovauth/ward.standaert@vliz.be/full_dat_outcomes")
# save(model_list, 
#      pr_predv, bg_predv, 
#      pr_pa, bg_pa, 
#      pr_coord, bg_coord, file = "full_dat_modelling_output.Rdata")



# library("usdm")
# plot(Variogram(st_list_NEA_cl[[1]][[1]][[1]], size=50))
# 
# moran(c(res),mat2listw(dists.inv), nrow(dists.inv), Szero(mat2listw(dists.inv)))
# MC <- moran.mc(c(res), mat2listw(dists.inv), nsim = 999, alternative="greater")
# plot(MC)
# 
# 
# 
# x_full <- full_dat %>% filter(scientificname == species$scientific[s])
# x_red <- reduced_dat %>% filter(scientificname == species$scientific[s])
# 
# tab_fun <- function(l) {
#   tab <- tibble(var = rep(names(l), 3),
#                 spat_temp = rep(c("Temporal_autocorrelation", "Temporal_autocorrelation", "Spatial_autocorrelation"), each = length(l)),
#                 mon_year = rep(c("month","year",NA), each = length(l))) %>%
#     unite(col = "id", var, spat_temp, mon_year, sep = ".", na.rm = TRUE, remove = FALSE)
#   dat <- unlist(l)
#   tab <- tab %>% 
#     rowwise() %>%
#     mutate(statistic = round(as.numeric(dat[which(str_detect(names(dat), id))][1]), 2),
#            p_value = round(as.numeric(dat[which(str_detect(names(dat), id))][3]), 2)) %>%
#     ungroup()
#   
#   return(tab)
# }
# 
# static_vars <- c("seabed_energy", "seabed_substrate","depth")
# dynamic_vars <- c("SST", "SSS", "windfarms", "Chl", "ZooPl", "max_SSV", "EuphD")
# SAC_list <- list()
# tmp_list <- list()
# 
# for (s in 1:nrow(species)) {
#   x_red <- reduced_dat %>% filter(scientificname == species$scientific[s]) %>%
#     mutate(x = lon,
#            y = lat,
#            year = Year,
#            month = Month,
#            day = Day,
#            seabed_energy = as.numeric(seabed_energy),
#            seabed_substrate = as.numeric(seabed_substrate)) %>%
#     select_if(function(col) length(unique(col)) > 1)
#     
#   x_full <- full_dat %>% filter(scientificname == species$scientific[s]) %>%
#     mutate(x = lon,
#            y = lat,
#            year = Year,
#            month = Month,
#            day = Day,
#            seabed_energy = as.numeric(seabed_energy),
#            seabed_substrate = as.numeric(seabed_substrate)) %>%
#     select_if(function(col) length(unique(col)) > 1)
#   
#   for (st in c("static", "dynamic")) {
#     v_list_st_d <- if(st == "static") static_vars else dynamic_vars
#     v_list_st_d_s <- v_list[[s]][which((v_list[[s]] %in% v_list_st_d) & (v_list[[s]] %in% names(x_red)))]
#     
#     cor_x_red <- spatiotemp_autocorr(occ.data = x_red, 
#                                      varname = v_list_st_d_s,
#                                      temporal.level = c("year","month"),
#                                      plot = FALSE)
#     cor_x_full <- spatiotemp_autocorr(occ.data = x_full, 
#                                       varname = v_list_st_d_s,
#                                       temporal.level = c("year","month"),
#                                       plot = FALSE)
#     a <- tab_fun(cor_x_full)
#     b <- tab_fun(cor_x_red)
#     
#     j <- left_join(a, b %>% dplyr::select(id, statistic, p_value), by = "id")
#     tmp_list[[st]] <- j
#   }
#   SAC_list[[s]] <- tmp_list
#   print(species$simple[s])
# }


# morans_I_list_spdep[[1]]   
# #After filtering:  Moran's I value of 0.118, expected -0.001; significant
# morans_I_list_spdep[[2]]   
# #After filtering:  Moran's I value of 0.0039, expected -0.0001; significant
# morans_I_list_spdep[[3]]   
# #After filtering:  Moran's I value of 0.095, expected -0.0065; significant
# morans_I_list_spdep[[4]]   
# #After filtering:  Moran's I value of 0.203, expected -0.007; significant

# Moran's I value of 0.073; expected -0.0004, variance 0.0001686635; sign

# 
# coordinates(cleandata_50) <- ~ Lat_UTM + Lon_UTM
# crs(cleandata_50) <- "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# 
# All_50_SPDF <- as(stack_50, "SpatialPixelsDataFrame")
# 
# pred_max <- predict(stack_50, maxnet_50, se.fit=TRUE, type="cloglog")
# 
# All_50_SPDF$GAM <- pred_gam$fit
# All_50_SPDF$GAM_se <- pred_gam$se.fit
# gam_pred <- stack(All_50_SPDF)
# 
# 
# pred_dat_max <- extract(pred_max, cleandata_50)
# 
# cleandata_50$pred_gam <- pred_dat_gam[,6]
# cleandata_50$pred_max <- pred_dat_max
# cleandata_50$pred_rf  <- pred_dat_rf
# cleandata_50$pred_ens <- pred_dat_ens
# 
# cleandata_50$res_gam <- cleandata_50$PresAbs - cleandata_50$pred_gam
# cleandata_50$res_max <- cleandata_50$PresAbs - cleandata_50$pred_max
# cleandata_50$res_rf  <- cleandata_50$PresAbs - cleandata_50$pred_rf
# cleandata_50$res_ens <- cleandata_50$PresAbs - cleandata_50$pred_ens
# df <- as.data.frame(cleandata_50)
# 
# ## folowing tutorial  on https://rstudio-pubs-static.s3.amazonaws.com/79757_3462f2bba8d745378545f0ea5bc38ee1.html
# 
# #calculate IDW
# dists <- as.matrix(dist(cbind(df$Lon_dec, df$Lat_dec)))
# dists.inv <- 1/dists
# diag(dists.inv) <- 0
# dists.inv[is.infinite(dists.inv)] <- 0 #remove infinite values
# 
# #Global Moran's I from ape package
# Moran.I(df$res_gam,dists.inv, scaled = TRUE, alternative = "greater")  #0.055, sign
# Moran.I(df$res_max,dists.inv, scaled = TRUE, alternative = "greater")  #0.341, sign
# Moran.I(df$res_rf,dists.inv, scaled = TRUE, alternative = "greater")   #-0.032, not sign
# 
# #Global Moran's I from spdep package are way slower (takes > 5 min for each)
# moran.test(df$res_gam,mat2listw(dists.inv)) #Moran's I value of 0.073; expected -0.0004, variance 0.0001686635; sign
# moran.test(df$res_max,mat2listw(dists.inv)) #Moran's I value of 0.778; expected -0.0004; sign
# moran.test(df$res_rf,mat2listw(dists.inv))  #Moran's I value of 0.002; expected -0.0004; not sign
# 
# #Local Moran's I based on k neirest neighbours
# coordinates(df) <- ~ Lon_UTM + Lat_UTM
# df2 <- remove.duplicates(df) #because error from knearneigh that there were duplicates --> to check
# w <- knn2nb(knearneigh(df2,k=4))
# moran.test(df2$res_gam,nb2listw(w)) #Moran's I value of 0.32
# moran.test(df2$res_max,nb2listw(w)) #Moran's I value of 0.48
# moran.test(df2$res_rf,nb2listw(w)) #Moran's I value of -0.02
# 
# 
# 
# library(gstat)
# corr.ran<-data.frame(res=df$res_max,lon=df$Lon_UTM,lat=df$Lat_UTM)
# #### check variogram
# coordinates(corr.ran)<-~lon+lat
# vario.ran<-variogram(res~lon+lat,data=corr.ran, cutoff = 1000) ;
# plot(vario.ran)

#6. variable importance ----

#make big tibble with all possible values per variable

load("~/INPUT/tib_all_v.rData")

# tib_all_v <- tibble(a = rep(1,3846528)) %>% dplyr::select()
# for (v in 1:length(v_list_all)) {
#   nm <- v_list_all[v]
#   tmp_x <- NULL
#   for (y in 1:length(2000:2020)) {
#     for (m in 1:12) {
#       tmp_st <-  st_list_NEA_cl[[y]][[m]]
#       tmp_x <- append(tmp_x, values(tmp_st[[which(str_detect(names(tmp_st), nm))]]))
#     }
#   }
#   tib_all_v <- tib_all_v %>% add_column(tmp_x)
#   colnames(tib_all_v)[v] <- v_list_all[v]
#   print(nm)
# }
# tib_all_v$Month <- rep(1:12, times = 21, each = ncell(st_list_NEA2[[1]][[1]]))
# tib_all_v$Year <- rep(2000:2020, times = 1, each = ncell(st_list_NEA2[[1]][[1]])*12)
# # tib_all_v$ac <- rep(0, nrow(tib_all_v))
# colnames(tib_all_v)[which(colnames(tib_all_v) %in% c("ene_char", "sub_char"))] <- c("seabed_energy", "seabed_substrate")
# 
# tib_all_v <- tib_all_v %>%
#   left_join(substr_lvl %>% mutate(seabed_substrate = as.numeric(seabed_substrate)), by = "seabed_substrate") %>%
#   left_join(energy_lvl %>% mutate(seabed_energy = as.numeric(seabed_energy)), by = "seabed_energy") %>%
#   left_join(month_lvl %>% mutate(Month = as.numeric(Month)), by = "Month") %>%
#   dplyr::select(-seabed_substrate, -seabed_energy, -Month)
# 
# setwd("/data/home/innovauth/ward.standaert@vliz.be/INPUT/")
# save(tib_all_v, file = "SAVE/tib_all_v.rData")

# tib_all_v <- tib_all_v %>%
#   left_join(substr_lvl %>% mutate(seabed_substrate = as.numeric(seabed_substrate)), by = "seabed_substrate") %>%
#   left_join(energy_lvl %>% mutate(seabed_energy = as.numeric(seabed_energy)), by = "seabed_energy") %>%
#   dplyr::select(-seabed_substrate, -seabed_energy)

tib_all_v <- tib_all_v %>%
  left_join(month_lvl, by = "mon_char") %>%
  mutate(Month_sin = sin(Month * (2*pi/12)),
         Month_cos = cos(Month * (2*pi/12)),
         quarter = ceiling(Month / 3),
         Quarter_sin = sin(quarter * (2*pi/4)),
         Quarter_cos = cos(quarter * (2*pi/4)))


cor_maxnet <- list()
for (s in 1:nrow(species)) {
  #1. get best model (based on AICc criterion)
  tic(paste0(species$simple[s], " - done"))
  
  #get model
  eval_res <- eval_res_list[[s]]
  opt.aicc <- eval.results(eval_res_list[[s]]) %>% filter(delta.AICc == 0)
  mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args[1])]]  
  
  #2. make original prediction
  df_pred <- rbind(pr_predv[[s]], bg_predv[[s]])
  if(all(c("seabed_energy", "seabed_substrate", "Month") %in% names(df_pred))) {
    df_pred <- df_pred %>%
      left_join(month_lvl, by = "Month") %>%
      left_join(energy_lvl, by = "seabed_energy") %>%
      left_join(substr_lvl, by = "seabed_substrate") %>%
      dplyr::select(-seabed_substrate, -seabed_energy, -Month)
  }
  else if(all(c("seabed_energy", "seabed_substrate") %in% names(df_pred))) {
    df_pred <- df_pred %>%
      left_join(energy_lvl, by = "seabed_energy") %>%
      left_join(substr_lvl, by = "seabed_substrate") %>%
      dplyr::select(-seabed_substrate, -seabed_energy)
  }
  else if (c("Month") %in% names(df_pred)) df_pred <- df_pred %>%
    left_join(month_lvl, by = "Month") %>%
    dplyr::select(-Month)
  
  
  
  prediction1 <- predict(mod, df_pred, se.fit=TRUE, type = "cloglog")
  
  #prepare cor matrix 
  cor_df_max <- data.frame(matrix(ncol = length(names(df_pred)), nrow= 50))
  colnames(cor_df_max) <- names(df_pred)
  
  for (v in 1:length(names(df_pred))) {
    #3. resample one variable
    var_name <- names(df_pred)[v]
    tic(var_name)
    for (i in 1:50) {
      df_pred_res <- df_pred
      ind_df <- which(names(df_pred_res) == var_name)
      ind_tib <- which(names(tib_all_v) == var_name)
      
      df_pred_res[,ind_df] <- sample(pull(na.omit(tib_all_v[,ind_tib])), nrow(df_pred_res))
      
      #sample each raster taking into account month & year - very slow!
      
      # df_pred_res <- df_pred
      # 
      # # Loop through unique combinations of Year and Month
      # unique_combinations <- df_pred_res %>%
      #   distinct(Year, Month)
      # 
      # for (i in 1:nrow(unique_combinations)) {
      #   year <- unique_combinations$Year[i]
      #   month <- unique_combinations$Month[i]
      #   
      #   # Get the rows that match the current Year and Month
      #   ind_rows <- which(df_pred_res$Month == month & df_pred_res$Year == year)
      #   
      #   # Filter and sample from tib_all_v
      #   sample_tib <- tib_all_v %>%
      #     filter(Year == year, Month == month) %>%
      #     dplyr::select({{var_name}}) %>%
      #     na.omit() %>%
      #     pull()
      #   
      #   n <- length(ind_rows)
      #   
      #   # Replace the values in df_pred_res
      #   df_pred_res[ind_rows, ind_df] <- sample(sample_tib, n)
      #   print(paste0(i,"/", nrow(unique_combinations)))
      # }
      
      #4. make new prediction
      prediction2 <- predict(mod, df_pred_res, se.fit=TRUE, type = "cloglog")
      cor_new_max <- cor(na.omit(prediction1), na.omit(prediction2))
      
      #calculate correlation between the two predictions
      cor_df_max[i, v] <- cor_new_max
    }
    toc()
  }
  #store in list
  cor_maxnet[[s]] <- sapply(cor_df_max, FUN = function(x) 1-mean(x, na.rm=T))
  
  toc()
}

save(cor_maxnet, file = "SAVE/NEA_variable_importance.Rdata")

# SAVE ----
load("SAVE/NEA_variable_importance.Rdata")

sort(cor_maxnet[[1]], decreasing = TRUE)
sort(cor_maxnet[[2]], decreasing = TRUE)
sort(cor_maxnet[[3]], decreasing = TRUE)
sort(cor_maxnet[[4]], decreasing = TRUE)


#7. plot variable importance & response curves ----
#plots using response.plot from maxnet package

for (s in 1:nrow(species)) {
  png(paste0("NEA_",species$simple[s], "_var_imp_response_curves.png"), width = 13, 
      height = 8, res = 100, units = 'in')
  
  #get model
  eval_res <- eval_res_list[[s]]
  opt.aicc <- eval.results(eval_res_list[[s]]) %>% filter(delta.AICc == 0)
  mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args[1])]]
  
  plot.new()
  #set plot layout
  layout.matrix <- matrix(c(0,1,0,2,3,4,5,6,7), nrow = 3, ncol = 3, byrow = TRUE)
  layout(mat = layout.matrix,
         heights = c(1.2,.8,.8),
         widths = c(1,1,1))
  
  #variable importance plot
  par(mar=c(5,6,4,2))
  barplot(sort(cor_maxnet[[s]][which(cor_maxnet[[s]] != 0)]), horiz = T, las = 1, xlab = "variable importance")
  
  #response curves
  v_imp <- sort(cor_maxnet[[s]], decreasing = TRUE)
  for (v in 1:length(which(v_imp > 0.05))) {
    v_name <- names(v_imp)[v]
    if(v_name == "ene_char") {
      par(mar=c(5,10,4,2))
      response.plot(mod, type = "cloglog", v = "ene_char", 
                    levels = na.omit(unlist(mod$levels$ene_char)[rev(match(energy_lvl$ene_char, unlist(mod$levels$ene_char)))]), las = 1, 
                    horiz = T, ylab= "", xlab = "Probability of occurrence - seabed energy",
                    cex.names = 0.8)
      next
    }
    if(v_name == "sub_char") {
      par(mar=c(5,14,4,2))
      response.plot(mod, type = "cloglog", v = "sub_char", 
                    levels = na.omit(unlist(mod$levels$sub_char)[rev(match(substr_lvl$sub_char, unlist(mod$levels$sub_char)))]), 
                    las = 1, horiz = T, ylab= "", xlab = "Probability of occurrence - seabed substrate",
                    cex.names = 0.8)
      next
    }
    if(v_name == "mon_char") {
      par(mar=c(5,10,4,2))
      response.plot(mod, type = "cloglog", v = "mon_char", 
                    levels = na.omit(unlist(mod$levels$mon_char)[rev(match(month_lvl$mon_char, unlist(mod$levels$mon_char)))]), 
                    las = 1, horiz = T, ylab= "", xlab = "Probability of occurrence - month",
                    cex.names = 0.8)
      next
    }
    min <- pr_predv[[s]] %>% dplyr::select(all_of(v_name)) %>% min
    max <- pr_predv[[s]] %>% dplyr::select(all_of(v_name)) %>% max
    
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    response.plot(mod, v_name, type = "cloglog", 
                  ylab = "Probability of occurrence",
                  min = min, max = max)
  }
  mtext(str_to_title(species$simple[s]), side = 3, line = - 2, outer = TRUE)
  mtext("Variable importance",           side = 3, line = - 4, outer = TRUE, cex = .8)
  mtext("Response curves",               side = 3, line = -29, outer = TRUE, cex = .8)
  dev.off()
}




#8. predictions NEA ----

## some tests to see how prediction works (esp. with categorical values) ----
# s <- 3
# 
# #get model
# eval_res <- eval_res_list[[s]]
# opt.aicc <- eval.results(eval_res_list[[s]]) %>% filter(delta.AICc == 0)
# mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]
# 
# names(mod$levels)
# plot_st_test <- stack(raster(matrix(c(sin(2*pi/6), sin(3*pi/6),
#                                       sin(4*pi/6), sin(5*pi/6)),
#                                     nrow=2, byrow = TRUE)),
#                       raster(matrix(c(cos(2*pi/6), cos(3*pi/6), 
#                                       cos(4*pi/6), cos(5*pi/6)),
#                                     nrow=2, byrow = TRUE)),
#                       raster(matrix(c(1,1,1,1),nrow=2)),
#                       raster(matrix(c(1,1,1,1),nrow=2)),
#                       raster(matrix(c(1,1,1,1),nrow=2)),
#                       raster(matrix(c(1,1,1,1),nrow=2)),
#                       raster(matrix(c(-100,-100,-100,-100),nrow=2)),
#                       raster(matrix(c(1,1,1,1),nrow=2)),
#                       raster(matrix(c(1,1,1,1),nrow=2)))
# names(plot_st_test) <- names(mod$levels)
# 
# substr_lvl
# mod$betas
# energy_lvl
# 
# pr_t <- predict(plot_st_test, mod, clamp=T, type="cloglog",
#                 factors = list(ene_char = factor(energy_lvl$seabed_energy,
#                                                  labels = energy_lvl$ene_char,
#                                                  levels = energy_lvl$seabed_energy),
#                                sub_char = factor(substr_lvl$seabed_substrate,
#                                                  labels = substr_lvl$sub_char,
#                                                  levels = substr_lvl$seabed_substrate)))
# plot(pr_t)

# plot IMDIS
# s <- 1
# 
# pr_geo <- stack()
# for (n in c(1,2)) {
#   m <- c(3,9)[n]
#   #get model
#   eval_res <- eval_res_list[[s]]
#   opt.aicc <- eval.results(eval_res_list[[s]]) %>% filter(delta.AICc == 0)
#   mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]
#   plot_st <- st_list_NEA_cl[[21]][[m]]
#   month_sin_r <- st_list_NEA_cl[[21]][[1]][[1]]
#   values(month_sin_r) <- sin(m * (2*pi/12))
#   names(month_sin_r) <- "Month_sin"
#   
#   month_cos_r <- st_list_NEA_cl[[21]][[1]][[1]]
#   values(month_cos_r) <- cos(m * (2*pi/12))
#   names(month_cos_r) <- "Month_cos"
#   
#   quarter_sin_r <- st_list_NEA_cl[[21]][[1]][[1]]
#   values(quarter_sin_r) <- sin(ceiling(m / 3) * (2*pi/12))
#   names(quarter_sin_r) <- "Quarter_sin"
#   
#   quarter_cos_r <- st_list_NEA_cl[[21]][[1]][[1]]
#   values(quarter_cos_r) <- cos(ceiling(m / 3) * (2*pi/12))
#   names(quarter_cos_r) <- "Quarter_cos"
#   
#   plot_st <- stack(plot_st, month_sin_r, month_cos_r, quarter_sin_r, quarter_cos_r)
#   
#   names(plot_st) <- str_remove_all(names(plot_st), "_\\d{4}_\\d{2}|_\\d{4}_\\d")
#   names(plot_st)[which(str_detect(names(plot_st),"seabed_energy"))] <- "ene_char"
#   names(plot_st)[which(str_detect(names(plot_st),"seabed_substrate"))] <- "sub_char"
#   
#   pr_geo <- stack(pr_geo, predict(plot_st, mod, clamp=T, type="cloglog", 
#                                   factors = list(ene_char = factor(energy_lvl$seabed_energy,
#                                                                    labels = energy_lvl$ene_char, 
#                                                                    levels = energy_lvl$seabed_energy),
#                                                  sub_char = factor(substr_lvl$seabed_substrate,
#                                                                    labels = substr_lvl$sub_char, 
#                                                                    levels = substr_lvl$seabed_substrate))))
#   names(pr_geo)[n] <- month.name[m]
# }
# 
# is_not_na <- !is.na(pr_geo[[1]])
# min(is_not_na)
# 
# mask_layer <- pr_geo[[1]]
# values(mask_layer) <- 1
# values(mask_layer)[which(is.na(values(pr_geo[[1]])))] <- 0
# 
# x_min <- min(coordinates(mask_layer)[which(values(mask_layer) == 1),1])
# x_max <- max(coordinates(mask_layer)[which(values(mask_layer) == 1),1])
# y_min <- min(coordinates(mask_layer)[which(values(mask_layer) == 1),2])
# y_max <- max(coordinates(mask_layer)[which(values(mask_layer) == 1),2])
# 
# pr_geo <- crop(pr_geo, extent(x_min,x_max,
#                               y_min,y_max))
# 
# levelplot(pr_geo, col.regions = colorRampPalette(rev(brewer.pal(10, 'RdYlBu'))), 
#           margin = FALSE, main = NULL, labels = FALSE, 
#           colorkey=list(title = expression(atop("Habitat","suitability")),
#                         cex.lab = 0.5, row=5, column=2, vjust=1.5, hjust= 0.8))


## prediction 2020 ----

for (y in 1:length(2000:2020)) {
  for (s in 1:nrow(species)) {
    pr_geo <- stack()
    for (m in 1:12) {
      #get model
      eval_res <- eval_res_list[[s]]
      opt.aicc <- eval.results(eval_res_list[[s]]) %>% filter(delta.AICc == 0)
      mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]
      plot_st <- st_list_NEA_cl[[y]][[m]]
      month_sin_r <- st_list_NEA_cl[[y]][[1]][[1]]
      values(month_sin_r) <- sin(m * (2*pi/12))
      names(month_sin_r) <- "Month_sin"
      
      month_cos_r <- st_list_NEA_cl[[y]][[1]][[1]]
      values(month_cos_r) <- cos(m * (2*pi/12))
      names(month_cos_r) <- "Month_cos"
      
      quarter_sin_r <- st_list_NEA_cl[[y]][[1]][[1]]
      values(quarter_sin_r) <- sin(ceiling(m / 3) * (2*pi/4))
      names(quarter_sin_r) <- "Quarter_sin"
      
      quarter_cos_r <- st_list_NEA_cl[[y]][[1]][[1]]
      values(quarter_cos_r) <- cos(ceiling(m / 3) * (2*pi/4))
      names(quarter_cos_r) <- "Quarter_cos"
      
      plot_st <- stack(plot_st, month_sin_r, month_cos_r, quarter_sin_r, quarter_cos_r)
      
      names(plot_st) <- str_remove_all(names(plot_st), "_\\d{4}_\\d{2}|_\\d{4}_\\d")
      names(plot_st)[which(str_detect(names(plot_st),"seabed_energy"))] <- "ene_char"
      names(plot_st)[which(str_detect(names(plot_st),"seabed_substrate"))] <- "sub_char"
      
      pr_geo <- stack(pr_geo, predict(plot_st, mod, clamp=T, type="cloglog", 
                                      factors = list(ene_char = factor(energy_lvl$seabed_energy,
                                                                       labels = energy_lvl$ene_char, 
                                                                       levels = energy_lvl$seabed_energy),
                                                     sub_char = factor(substr_lvl$seabed_substrate,
                                                                       labels = substr_lvl$sub_char, 
                                                                       levels = substr_lvl$seabed_substrate))))
      names(pr_geo)[m] <- month.abb[m]
    }
    
    is_not_na <- !is.na(pr_geo[[1]])
    min(is_not_na)
    
    mask_layer <- pr_geo[[1]]
    values(mask_layer) <- 1
    values(mask_layer)[which(is.na(values(pr_geo[[1]])))] <- 0
    
    x_min <- min(coordinates(mask_layer)[which(values(mask_layer) == 1),1])
    x_max <- max(coordinates(mask_layer)[which(values(mask_layer) == 1),1])
    y_min <- min(coordinates(mask_layer)[which(values(mask_layer) == 1),2])
    y_max <- max(coordinates(mask_layer)[which(values(mask_layer) == 1),2])
    
    pr_geo <- crop(pr_geo, extent(x_min,x_max,
                                  y_min,y_max))
    

    #save picture
    png(paste0("output_rasters_NEA/NEA_",c(2000:2020)[y],"_",species$simple[s], "_prediction.png"), width = 8, height = 9, res = 100, units = 'in')
    print(levelplot(pr_geo, col.regions = colorRampPalette(rev(brewer.pal(10, 'RdYlBu'))), 
                    margin = FALSE, main = paste0("HSI ", species$scientific[s])))
    dev.off()    
    
    #save rasters
    
    if(y == 21) {
      for(i in 1:nlayers(pr_geo)) {
        writeRaster(pr_geo[[i]], paste0("output_rasters_NEA/NEA_",species$simple[s],"_adult_",c(2000:2020)[y],"_",i,".tif"))
        }
    }
    
    print(paste0(c(2000:2020)[y],"_",species$simple[s]))
  }
}


#9. HSI evolution over years ----
for (m in 1:12) {
  pr_geo_y <- stack()
  for (y in 1:length(2000:2020)) {
    #get model
    eval_res <- eval_res_list[[1]]
    opt.aicc <- eval.results(eval_res_list[[1]]) %>% filter(delta.AICc == 0)
    mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]
    plot_st <- st_list_NEA_cl[[y]][[m]]
    month_sin_r <- st_list_NEA_cl[[y]][[1]][[1]]
    values(month_sin_r) <- sin(m * (2*pi/12))
    names(month_sin_r) <- "Month_sin"
    
    month_cos_r <- st_list_NEA_cl[[y]][[1]][[1]]
    values(month_cos_r) <- cos(m * (2*pi/12))
    names(month_cos_r) <- "Month_cos"
    
    quarter_sin_r <- st_list_NEA_cl[[y]][[1]][[1]]
    values(quarter_sin_r) <- sin(ceiling(m / 3) * (2*pi/4))
    names(quarter_sin_r) <- "Quarter_sin"
    
    quarter_cos_r <- st_list_NEA_cl[[y]][[1]][[1]]
    values(quarter_cos_r) <- cos(ceiling(m / 3) * (2*pi/4))
    names(quarter_cos_r) <- "Quarter_cos"
    
    plot_st <- stack(plot_st, month_sin_r, month_cos_r, quarter_sin_r, quarter_cos_r)
    
    names(plot_st) <- str_remove_all(names(plot_st), "_\\d{4}_\\d{2}|_\\d{4}_\\d")
    names(plot_st)[which(str_detect(names(plot_st),"seabed_energy"))] <- "ene_char"
    names(plot_st)[which(str_detect(names(plot_st),"seabed_substrate"))] <- "sub_char"
    
    pr_geo_y <- stack(pr_geo_y, predict(plot_st, mod, clamp=T, type="cloglog", 
                                        factors = list(ene_char = factor(energy_lvl$seabed_energy,
                                                                         labels = energy_lvl$ene_char, 
                                                                         levels = energy_lvl$seabed_energy),
                                                       sub_char = factor(substr_lvl$seabed_substrate,
                                                                         labels = substr_lvl$sub_char, 
                                                                         levels = substr_lvl$seabed_substrate))))
    names(pr_geo_y)[y] <- c(2000:2020)[y]
  }
  
  # average & coefficient of variation NEA
  av_NEA <- stackApply(pr_geo_y, indices =  rep(1,nlayers(pr_geo_y)), mean, na.rm = T)
  sd_NEA <- stackApply(pr_geo_y, indices =  rep(1,nlayers(pr_geo_y)), sd, na.rm = T)
  uncert_NEA <- sd_NEA/av_NEA*100
  
  writeRaster(av_NEA, paste0("output_rasters_NEA/NEA_herring_adult_average_HSI_",m,".tif"), overwrite = T)
  writeRaster(sd_NEA, paste0("output_rasters_NEA/NEA_herring_adult_sd_HSI_",m,".tif"), overwrite = T)
  writeRaster(uncert_NEA, paste0("output_rasters_NEA/NEA_herring_adult_CV_HSI_",m,".tif"), overwrite = T)
  
  # variability BPNS
  pr_geo_y <- crop(pr_geo_y, bpns_shp)
  pr_geo_y <- mask(pr_geo_y, bpns_shp)
  
  df <- values(pr_geo_y) %>%
    as.data.frame() %>%
    na.omit() %>%
    gather(key = year, value = "value") %>%
    mutate(year = str_remove(year,"X"))
  
  
  ggplot(df) +
    geom_boxplot(aes(year, value)) +
    scale_y_continuous(limits = c(0,1)) +
    xlab("Year") + 
    ylab("Habitat suitability index") +
    theme_bw() +
    theme(legend.title = element_blank(), 
          axis.text.y = element_text(size=10, face = "plain", colour = "black"),
          axis.text.x = element_text(size=10, face = "plain", colour = "black"),
          axis.title.x = element_text(size=10, face = "bold", colour = "black"),
          axis.title.y = element_text(size=10, face = "bold", colour = "black"),
          legend.text = element_text(size=10, face = "bold", colour = "black"))
  
  ggsave(paste0("evolution/adult_yearly_evolution_HSI_", m,".png"), width = 10, height = 4)
  
  print(m)
}


#10. HSI evolution over months ----
df_out <- data.frame(month = rep(month.abb, each = 9))
for (y in 1:length(2000:2020)) {
  pr_geo_y <- stack()
  for (m in 1:12) {
    #get model
    eval_res <- eval_res_list[[1]]
    opt.aicc <- eval.results(eval_res_list[[1]]) %>% filter(delta.AICc == 0)
    mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]
    plot_st <- st_list_NEA_cl[[y]][[m]]
    month_sin_r <- st_list_NEA_cl[[y]][[1]][[1]]
    values(month_sin_r) <- sin(m * (2*pi/12))
    names(month_sin_r) <- "Month_sin"
    
    month_cos_r <- st_list_NEA_cl[[y]][[1]][[1]]
    values(month_cos_r) <- cos(m * (2*pi/12))
    names(month_cos_r) <- "Month_cos"
    
    quarter_sin_r <- st_list_NEA_cl[[y]][[1]][[1]]
    values(quarter_sin_r) <- sin(ceiling(m / 3) * (2*pi/4))
    names(quarter_sin_r) <- "Quarter_sin"
    
    quarter_cos_r <- st_list_NEA_cl[[y]][[1]][[1]]
    values(quarter_cos_r) <- cos(ceiling(m / 3) * (2*pi/4))
    names(quarter_cos_r) <- "Quarter_cos"
    
    plot_st <- stack(plot_st, month_sin_r, month_cos_r, quarter_sin_r, quarter_cos_r)
    
    names(plot_st) <- str_remove_all(names(plot_st), "_\\d{4}_\\d{2}|_\\d{4}_\\d")
    names(plot_st)[which(str_detect(names(plot_st),"seabed_energy"))] <- "ene_char"
    names(plot_st)[which(str_detect(names(plot_st),"seabed_substrate"))] <- "sub_char"
    
    pr_geo_y <- stack(pr_geo_y, predict(plot_st, mod, clamp=T, type="cloglog", 
                                        factors = list(ene_char = factor(energy_lvl$seabed_energy,
                                                                         labels = energy_lvl$ene_char, 
                                                                         levels = energy_lvl$seabed_energy),
                                                       sub_char = factor(substr_lvl$seabed_substrate,
                                                                         labels = substr_lvl$sub_char, 
                                                                         levels = substr_lvl$seabed_substrate))))
    names(pr_geo_y)[m] <- month.abb[m]
  }
  
  pr_geo_y <- crop(pr_geo_y, bpns_shp)
  pr_geo_y <- mask(pr_geo_y, bpns_shp)
  
  df <- values(pr_geo_y) %>%
    as.data.frame() %>%
    na.omit() %>%
    gather(key = month, value = "value") %>%
    mutate(month = str_remove(month,"X"))
  
  df_out <- cbind(df_out, df[,2])
  names(df_out)[y+1] <- as.character(c(2000:2020)[y])
  
  ggplot(df) +
    geom_boxplot(aes(month, value)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(limits = month.abb) +
    xlab("Month") + 
    ylab("Habitat suitability index") +
    theme_bw() +
    theme(legend.title = element_blank(), 
          axis.text.y = element_text(size=10, face = "plain", colour = "black"),
          axis.text.x = element_text(size=10, face = "plain", colour = "black"),
          axis.title.x = element_text(size=10, face = "bold", colour = "black"),
          axis.title.y = element_text(size=10, face = "bold", colour = "black"),
          legend.text = element_text(size=10, face = "bold", colour = "black"))
  
  ggsave(paste0("evolution/adult_monthly_evolution_HSI_", c(2000:2020)[y],".png"), width = 10, height = 4)
  
  print(c(2000:2020)[y])
}

write.csv(df_out, "evolution/evolution_BPNS_values.csv")

#11. average monthly HSI ----
df_evol <- read.csv("evolution/evolution_BPNS_values.csv")

df_evol <- df_evol %>%
  select(-X) %>%
  pivot_longer(cols = c(-month))

ggplot(df_evol) +
  geom_boxplot(aes(month, value)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(limits = month.abb) +
  xlab("Month") + 
  ylab("Habitat suitability index") +
  theme_bw() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size=10, face = "plain", colour = "black"),
        axis.text.x = element_text(size=10, face = "plain", colour = "black"),
        axis.title.x = element_text(size=10, face = "bold", colour = "black"),
        axis.title.y = element_text(size=10, face = "bold", colour = "black"),
        legend.text = element_text(size=10, face = "bold", colour = "black"))

ggsave(paste0("evolution/adult_monthly_evolution_HSI_average.png"), width = 10, height = 4)

#12. average yearly HSI ----
df_evol <- read.csv("evolution/evolution_BPNS_values.csv")

df_evol <- df_evol %>%
  select(-X, -month) %>%
  gather(key = "year", value = "value") 

df_evol$year <- str_remove(df_evol$year, "X")

ggplot(df_evol) +
  geom_boxplot(aes(year, value)) +
  scale_y_continuous(limits = c(0,1)) +
  xlab("Year") + 
  ylab("Habitat suitability index") +
  theme_bw() +
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(size=10, face = "plain", colour = "black"),
        axis.text.x = element_text(size=10, face = "plain", colour = "black"),
        axis.title.x = element_text(size=10, face = "bold", colour = "black"),
        axis.title.y = element_text(size=10, face = "bold", colour = "black"),
        legend.text = element_text(size=10, face = "bold", colour = "black"))

ggsave(paste0("evolution/adult_yearly_evolution_HSI_average.png"), width = 10, height = 4)

#13. HSI over latitude ----
lat_v_all <- unique(coordinates(st_list_NEA_cl[[1]][[1]])[,2])
df <- data.frame(lat = rep(lat_v_all,21),
                 val = rep(NA),
                 year = c(rep(2000:2020, each = length(lat_v_all))))

for (y in 1:length(2000:2020)) {
  m <- 2 #for january
  
  eval_res <- eval_res_list[[1]]
  opt.aicc <- eval.results(eval_res_list[[1]]) %>% filter(delta.AICc == 0)
  mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]
  plot_st <- st_list_NEA_cl[[y]][[m]]
  month_sin_r <- st_list_NEA_cl[[y]][[1]][[1]]
  values(month_sin_r) <- sin(m * (2*pi/12))
  names(month_sin_r) <- "Month_sin"
  
  month_cos_r <- st_list_NEA_cl[[y]][[1]][[1]]
  values(month_cos_r) <- cos(m * (2*pi/12))
  names(month_cos_r) <- "Month_cos"
  
  quarter_sin_r <- st_list_NEA_cl[[y]][[1]][[1]]
  values(quarter_sin_r) <- sin(ceiling(m / 3) * (2*pi/4))
  names(quarter_sin_r) <- "Quarter_sin"
  
  quarter_cos_r <- st_list_NEA_cl[[y]][[1]][[1]]
  values(quarter_cos_r) <- cos(ceiling(m / 3) * (2*pi/4))
  names(quarter_cos_r) <- "Quarter_cos"
  
  plot_st <- stack(plot_st, month_sin_r, month_cos_r, quarter_sin_r, quarter_cos_r)
  
  names(plot_st) <- str_remove_all(names(plot_st), "_\\d{4}_\\d{2}|_\\d{4}_\\d")
  names(plot_st)[which(str_detect(names(plot_st),"seabed_energy"))] <- "ene_char"
  names(plot_st)[which(str_detect(names(plot_st),"seabed_substrate"))] <- "sub_char"
  
  pr_geo <- predict(plot_st, mod, clamp=T, type="cloglog", 
                    factors = list(ene_char = factor(energy_lvl$seabed_energy,
                                                     labels = energy_lvl$ene_char, 
                                                     levels = energy_lvl$seabed_energy),
                                   sub_char = factor(substr_lvl$seabed_substrate,
                                                     labels = substr_lvl$sub_char, 
                                                     levels = substr_lvl$seabed_substrate)))
  
  for (l in 1:length(lat_v_all)) {
    v <- mean(values(pr_geo)[which(coordinates(pr_geo)[,2] == lat_v_all[l])], na.rm = T)
    df$val[l + (y-1)*length(lat_v_all)] <- v
  }
  print(y)
}

ggplot(df, aes(x = lat, y = val, group = year)) +
  geom_line(aes(colour = year)) +
  scale_x_continuous(limits = c(52,60)) +
  scale_color_continuous(type = "viridis") +
  scale_y_continuous(limits = c(0.4,0.6)) +
  labs(x = "Latitude",
       y = "Suitability index") +
  theme_minimal()
