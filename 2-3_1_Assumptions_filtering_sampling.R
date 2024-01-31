pckgs <- c("raster","sp","proj4","ncdf4","car","mgcv","dismo","rJava","ENMeval","randomForest", "VSURF",
           "randomForestSRC","ggRandomForests","ggplot2","boot","gstat","mgcv","SDMtune","ENMeval",
           "tidyverse","mapdata","base","tidync", "sf", "usdm", "mapview", "tictoc", "ape", "spdep",
           "spThin", "StatMatch")
lapply(pckgs, require, character.only = TRUE)
rm(pckgs)
setwd("DATA")

#1. Read & organize data -----

##1.1. Biological data ----
###1.1.1. Adult data ----
load("2.PREPROCESSED/biology/her_mack_tws_seab_adults.Rdata")

#read spawning data
herring_sp <- read.csv("2.PREPROCESSED/biology/OBIS_larvae_herring_2000-2020.csv")
mackerel_sp <- read.csv("2.PREPROCESSED/biology/OBIS_larvae_mackerel_2000-2020.csv")

dati_sp <- rbind(herring_sp %>% dplyr::select(lon, lat, Year = year, Month = month, Day = day) %>%
                   mutate(scientificname = "Clupea harengus"),
                 mackerel_sp %>% dplyr::select(lon, lat, Year = year, Month = month, Day = day) %>%
                   mutate(scientificname = "Scomber scombrus"))
rm(herring_sp, mackerel_sp)


#### Remove spatial outliers ----
#1. compiled plot
x <- data.frame(lon = dati_ad$lon, 
                lat = dati_ad$lat, 
                species = dati_ad$scientificname)
coordinates(x) <- ~lon+lat
crs(x) <- "epsg:4326"
mapview(x, zcol = "species")
#remove Baltic Sea observations; at longitude 9.3656°, according to marine regions
table(dati_ad$scientificname)
dati_ad <- dati_ad %>%
  filter(lon <= 9.3656)
table(dati_ad$scientificname)

#remove outlying observations of mackerel in northeast
dati_ad <- dati_ad %>%
  filter(!(lon < -12 & lat > 55))
table(dati_ad$scientificname)

#2. per species
mapview(dati_ad %>% filter(scientificname == "Clupea harengus") %>%
          select(lon) %>% pull(), 
        dati_ad %>% filter(scientificname == "Clupea harengus") %>%
          select(lat) %>% pull(), 
        crs = "epsg:4326")
#remove points in southwest (lat < 48) & (lat < 50 & lon < -10)
mapview(dati_ad %>% filter(scientificname == "Scomber scombrus") %>%
          select(lon) %>% pull(), 
        dati_ad %>% filter(scientificname == "Scomber scombrus") %>%
          select(lat) %>% pull(), 
        crs = "epsg:4326")
mapview(dati_ad %>% filter(scientificname == "Alosa fallax") %>%
          select(lon) %>% pull(), 
        dati_ad %>% filter(scientificname == "Alosa fallax") %>%
          select(lat) %>% pull(), 
        crs = "epsg:4326")
#remove one point west of Ireland (lon < -10)
mapview(dati_ad %>% filter(scientificname == "Dicentrarchus labrax") %>%
          select(lon) %>% pull(), 
        dati_ad %>% filter(scientificname == "Dicentrarchus labrax") %>%
          select(lat) %>% pull(), 
        crs = "epsg:4326")
#remove 2 solitary points off of Spain (lon < -5 & lat < 45)
#remove 3 solitary points in east (lon > 5)

table(dati_ad$scientificname)
dati_ad <- dati_ad %>% 
  filter(!(lat < 48 & scientificname == "Clupea harengus")) %>%
  filter(!(lat < 50 & lon < -10 & scientificname == "Clupea harengus")) %>%
  filter(!(lon < -10 & scientificname == "Alosa fallax")) %>%
  filter(!(lon < -5 & lat < 45 & scientificname == "Dicentrarchus labrax")) %>%
  filter(!(lon > 5 & scientificname == "Dicentrarchus labrax"))
table(dati_ad$scientificname)


#### Remove temporal outliers ----
#first remove months that are not represented in surveys
table(dati_ad$scientificname, dati_ad$Month)
#                         1    2    3    4    6    7    8    9   10   11   12
# Alosa fallax          108  266   11    0    0   18   99   50  114   27    2
# Clupea harengus      1542 3493 1472   51    2  674 2685  689  887 1346  271
# Dicentrarchus labrax   72   70   59    9    0   79   56  234  834  252   23
# Scomber scombrus      243  522  427   10    3  572 2637 1058 2018 1462  238

# remove quarter 2, months 4-6 and months 3 and 12 for Alosa fallax
dati_ad <- dati_ad %>% 
  filter(!(Month %in% c(4,6))) %>%
  filter(!(scientificname == "Alosa fallax" & Month %in% c(3,12)))

table(dati_ad$scientificname, dati_ad$Month)
#                         1    2    3    7    8    9   10   11   12
# Alosa fallax          108  266    0   18   99   50  114   27    0
# Clupea harengus      1542 3493 1472  674 2685  689  887 1346  271
# Dicentrarchus labrax   72   70   59   79   56  234  834  252   23
# Scomber scombrus      243  522  427  572 2637 1058 2018 1462  238 


###1.1.2. Larval data ----
# as a proxy for spawning areas, abbreviation sp

#### Remove spatial outliers ----
mapview(dati_sp %>% filter(scientificname == "Clupea harengus") %>% select(lon) %>% pull,
        dati_sp %>% filter(scientificname == "Clupea harengus") %>% select(lat) %>% pull,
        crs = "epsg:4326")
mapview(dati_sp %>% filter(scientificname == "Scomber scombrus") %>% select(lon) %>% pull,
        dati_sp %>% filter(scientificname == "Scomber scombrus") %>% select(lat) %>% pull,
        crs = "epsg:4326")
#no outliers

#### Remove temporal outliers ----
#remove months that are underrepresented in surveys (month 8 for herring)
table(dati_sp$scientificname, dati_sp$Month)
#                      1     2     3     4     5     6     7     8     9    10    12
# Clupea harengus   2293     0     0     0     0     0     0     0  4480    99  1030
# Scomber scombrus    70   765  6423 12325 13531 11426  1776     1     0     0     0

dati_sp <- dati_sp %>% filter(!(Month %in% c(8)))

table(dati_sp$scientificname, dati_sp$Month)
#                      1     2     3     4     5     6     7     9    10    12
# Clupea harengus   2293     0     0     0     0     0     0  4480    99  1030
# Scomber scombrus    70   765  6423 12325 13531 11426  1776     0     0     0


#### Add ICES Area 27 to larvae occurrences ----
# ICES areas derived from https://gis.ices.dk on 18/12/2023
ices_shp <- st_read("1.DOWNLOAD/ICES_areas/ICES_Areas_20160601_cut_dense_3857.shp")
mapview(ices_shp)

ices_shp_t <- st_transform(ices_shp, "epsg:4326") %>%
  st_make_valid() %>%
  st_crop(c(xmin = -26, xmax = 40, ymin = 30, ymax = 70)) %>%
  st_make_valid()

coordinates_sf_sp <- st_as_sf(data.frame(lon = dati_sp$lon, lat = dati_sp$lat), coords = c("lon", "lat"), crs = st_crs(ices_shp_t))

tmp_ind <- apply(st_intersects(coordinates_sf_sp, ices_shp_t, sparse = FALSE), 1, function(x) any(x))

f_t <- function(df) {
  out <- which(df)
  if(length(out) == 0) NA
  else ices_shp_t$Area_27[out]
}

dati_sp$Area_27 <- apply(st_intersects(coordinates_sf_sp, ices_shp_t, sparse = FALSE), 1, f_t)
head(dati_sp)

##1.2. Environmental variables ----
# only run once
### create stack list from preprocessed files ----
## results in stacked list: st_list_NEA[[y]][[m]][[v]] with y for year, m for month and v for variable
# fls_NEA <- tibble(full_name = list.files("2.PREPROCESSED/environmental_variables", pattern = ".grd$"),
#                   shrt_name = str_remove_all(full_name, ".grd$"))
# 
# st_list_NEA <- list()
# tmp_st_list <- list()
# for (y in 1:length(2000:2020)) {
#   for (m in 1:12) {
#     tmp_st_m <- stack()
#     for (v in 1:nrow(fls_NEA)) {
#       tmp_st <- stack(paste("2.PREPROCESSED/environmental_variables/", fls_NEA$full_name[v], sep = ""))
#       if (nlayers(tmp_st) == 1) tmp_r <- tmp_st
#       else {
#           tmp_nms <- tibble(name = names(tmp_st),
#                             year = str_extract(name, "\\d{4}"),
#                             month = str_remove(str_extract(name, "_\\d{2}|_\\d"), "_"))
#           tmp_ind <- which(tmp_nms$month == m & tmp_nms$year == as.character(c(2000:2020))[y])
#           tmp_r <- tmp_st[[tmp_ind]]
#           }
#       names(tmp_r) <- paste(fls_NEA$shrt_name[v], c(2000:2020)[y], m, sep = "_")
#       tmp_st_m <- stack(tmp_st_m, tmp_r)
#     }
#     tmp_st_list[[m]] <- tmp_st_m
#   }
#   st_list_NEA[[y]] <- tmp_st_list
#   print(paste("processed", c(2000:2020)[y]))
# }
# 
# rm(tmp_st, tmp_st_m, tmp_st_list, tmp_nms, tmp_ind, tmp_r)
# 
# st_list_NEA[[1]][[1]][[1]]
# st_list_NEA[[1]][[1]][[2]]
# st_list_NEA[[1]][[1]][[1]]
# st_list_NEA[[1]][[2]][[1]]
# st_list_NEA[[1]][[1]][[1]]
# st_list_NEA[[2]][[1]][[1]]
# 
# save(st_list_NEA, file = "SAVE/stacked_vars_not_cropped.Rdata")

# load files directly
load("SAVE/stacked_vars_not_cropped.Rdata")

# only run once
### create an average and standard deviation of layers over 2000-2020 for each environmental variable ----
# for (v in 1:nrow(fls_NEA)) {
#   tmp_st <- stack(paste("2.PREPROCESSED/environmental_variables/", fls_NEA$full_name[v], sep = ""))
#   # no need to do this for static variables
#   if (nlayers(tmp_st) == 1) next
#   else {
#     tmp_nms <- tibble(name = names(tmp_st),
#                       year = str_extract(name, "\\d{4}"),
#                       month = str_remove(str_extract(name, "_\\d{2}|_\\d"), "_"))
#     for (m in 1:12) {
#       tmp_ind <- which(tmp_nms$month == m)
#       av_NEA <- stackApply(tmp_st[[tmp_ind]], indices =  rep(1, length(tmp_ind)), mean, na.rm = T)
#       sd_NEA <- stackApply(tmp_st[[tmp_ind]], indices =  rep(1, length(tmp_ind)), sd, na.rm = T)
#       writeRaster(av_NEA,
#                   paste0("2.PREPROCESSED/environmental_variables/average_standard_deviation/NEA_average_",
#                          fls_NEA$shrt_name[v],"_",
#                          m,".tif"), overwrite = T)
#       writeRaster(sd_NEA,
#                   paste0("2.PREPROCESSED/environmental_variables/average_standard_deviation/NEA_sd_",
#                          fls_NEA$shrt_name[v],"_",
#                          m,".tif"), overwrite = T)
#     }
#   }
# }

### crop of layers ----
# first general crop of layers to their common extent
plot(st_list_NEA[[1]][[1]])

## crop raster towards the ICES regions where observations are present
## also crop to common non-NA values of all rasters (common mask)

common_masks <- list()
for (y in 1:length(2000:2020)) {
  masks <- lapply(st_list_NEA[[y]], function(month_layer) {
    is_not_na <- !is.na(month_layer)
    min(is_not_na)
  })
  
  # Reduce the list of masks to get a common mask for the current year
  common_mask <- Reduce(`&`, masks)
  
  # Add the common mask to the list of common masks
  common_masks[[y]] <- common_mask
  print(y)
}

final_common_mask <- Reduce(`&`, common_masks)
values(final_common_mask)[which(values(final_common_mask) == 0)] <- NA

st_list_NEA2 <- list()
tmp_st_list <- stack()
for (y in 1:length(2000:2020)) {
  tmp_st_list <- lapply(st_list_NEA[[y]], mask, mask = final_common_mask)
  st_list_NEA2[[y]] <- tmp_st_list
  print(y)
}
rm(tmp_st_list)

plot(st_list_NEA2[[1]][[1]])

all_vars <- str_remove_all(names(st_list_NEA2[[1]][[1]]) , "_\\d{4}_\\d{2}|_\\d{4}_\\d")

### check for multi-collinearity ----
#### per month ----
# VIF calculating function
f_vif <- function(stack_list, retained_vars) {
  vif_tab <- data.frame(Var = retained_vars)
  for (y in 1:length(2000:2020)) {
    for (m in 1:12) {
      tmp_st <- st_list_NEA2[[y]][[m]][[which(str_detect(names(st_list_NEA2[[y]][[m]]), paste(retained_vars, collapse = "|")))]]
      st_df <- as.data.frame(as(tmp_st, "SpatialPixelsDataFrame")) %>%
        select(-x,-y)
      
      tmp_vif <- vif(st_df)
      tmp_vif_tab <- tmp_vif %>%
        mutate(Var = str_remove_all(Variables , "_\\d{4}_\\d{2}|_\\d{4}_\\d")) %>%
        select(Var, VIF)
      colnames(tmp_vif_tab) <- c("Var", paste("VIF",c(2000:2020)[y],m,sep = "_"))
      vif_tab <- left_join(vif_tab, tmp_vif_tab, by = "Var")
    }
    print(y)
  }
  vif_tab$av  <- rowMeans(vif_tab[,-1], na.rm = T)
  vif_tab$max <- apply(vif_tab[,-1], 1, max, na.rm = TRUE)
  return(vif_tab %>%
           relocate(Var, av, max) %>%
           arrange(desc(max)))
}


#including all variables
vif_t1 <- f_vif(st_list_NEA2, all_vars)
vif_t1$Var[apply(vif_t1[,-c(1,2)] > 10, 1, any)]
vif_t1[,1:3]
# chlorophyll and phytoplankton are highly correlated --> remove chlorophyll

retained_vars <- all_vars[-which(all_vars == "Chl")]
vif_t2 <- f_vif(st_list_NEA2, retained_vars)
any(vif_t2[,-c(1,2)] > 10)
vif_t2$Var[apply(vif_t2[,-c(1,2)] > 10, 1, any)]
vif_t2[,1:3]
# temperature and dissolved oxygen are highly correlated --> remove dissolved oxygen

retained_vars <- all_vars[-which(all_vars == "Chl" | all_vars == "DO" )]
vif_t3 <- f_vif(st_list_NEA2, retained_vars)
vif_t3[,1:3]
# euphotic depth has high max VIF --> remove euphotic depth

retained_vars <- all_vars[-which(all_vars == "Chl" | all_vars == "DO" | all_vars == "EuphD")]
vif_t4 <- f_vif(st_list_NEA2, retained_vars)
vif_t4[,1:3]
#conclusion: remove Chl, DO and EuphD

#remove these layers from stack
v_list <- all_vars[-which(all_vars == "Chl" | all_vars == "DO" | all_vars == "EuphD")]

st_list_NEA_cl <- st_list_NEA2
for (y in 1:length(2000:2020)) {
  for (m in 1:12) {
    st_list_NEA_cl[[y]][[m]] <- st_list_NEA2[[y]][[m]][[which(str_detect(names(st_list_NEA2[[y]][[m]]), paste(v_list, collapse = "|")))]]
  }
}

#2. Sample environmental variables with biological data ----

##2.1 Adult data ----
### remove duplicates in space & time ----
dati_ad <- dati_ad %>% 
  distinct(Year, Month, lon, lat, scientificname, .keep_all = TRUE)

### sampling ----
coords_dat <- dati_ad %>% select(lon, lat, Year, Month, Day, scientificname, Area_27)
full_dat <- data.frame()
for (y in 1:length(2000:2020)) {
  for (m in 1:12) {
    tmp_coords_dat <- coords_dat %>%
      filter(Year == c(2000:2020)[y],
             Month == m)
    tmp_extract <- raster::extract(st_list_NEA_cl[[y]][[m]], 
                                   as.data.frame(tmp_coords_dat[,c(1,2)]), method = "simple")
    colnames(tmp_extract) <- colnames(tmp_extract) %>%
      str_remove_all("_\\d{4}_\\d{2}|_\\d{4}_\\d")
    tmp_df <- data.frame(tmp_coords_dat, tmp_extract)
    full_dat <- rbind(full_dat, tmp_df)
  }
  print(y)
}
rm(tmp_extract, tmp_coords_dat, tmp_df)

# investigate observations without corresponding environmental values (NA-values) ----
any(is.na(full_dat$Gear)) #no NA values
any(is.na(full_dat$Ship)) #no NA values
na_dat <- full_dat %>% filter(is.na(depth) | is.na(max_SSV) | is.na(Phyto) | 
                                is.na(seabed_energy) | is.na(seabed_substrate) | is.na(SSS) | 
                                is.na(SST) | is.na(windfarms) | is.na(ZooPl))

head(na_dat)
mapview::mapview(na_dat$lon, na_dat$lat, crs = CRS("epsg:4326"), legend = FALSE) +
  mapview::mapview(st_list_NEA_cl[[1]][[3]][[1]], legend = FALSE)

na_dat2 <- na_dat %>% drop_na(Phyto)
head(na_dat2) ##all variables give NA-values at same points (= check that cropping was done well)

rm(na_dat, na_dat2)

#remove observations without corresponding environmental values
full_dat <- full_dat %>%
  drop_na()

##2.2 Larval data ----
### remove duplicates ----
dati_sp <- dati_sp %>% 
  distinct(Year, Month, lon, lat, scientificname, .keep_all = TRUE)

### sampling ----
coords_dat_sp <- dati_sp %>% select(lon, lat, Year, Month, Day, scientificname, Area_27)
full_dat_sp <- data.frame()
for (y in 1:length(2000:2020)) {
  for (m in 1:12) {
    tmp_coords_dat <- coords_dat_sp %>%
      filter(Year == c(2000:2020)[y],
             Month == m)
    tmp_extract <- raster::extract(st_list_NEA_cl[[y]][[m]], 
                                   as.data.frame(tmp_coords_dat[,c(1,2)]), method = "simple")
    colnames(tmp_extract) <- colnames(tmp_extract) %>%
      str_remove_all("_\\d{4}_\\d{2}|_\\d{4}_\\d")
    tmp_df <- data.frame(tmp_coords_dat, tmp_extract)
    full_dat_sp <- rbind(full_dat_sp, tmp_df)
  }
  print(y)
}
rm(tmp_extract, tmp_coords_dat, tmp_df)

#investigate observations without corresponding environmental values (NA-values)
na_dat <- full_dat_sp %>% filter(is.na(depth) | is.na(max_SSV) | is.na(Phyto) | 
                                is.na(seabed_energy) | is.na(seabed_substrate) | is.na(SSS) | 
                                is.na(SST) | is.na(windfarms) | is.na(ZooPl))

head(na_dat)
mapview::mapview(na_dat$lon, na_dat$lat, crs = CRS("epsg:4326"), legend = FALSE) +
  mapview::mapview(st_list_NEA_cl[[1]][[3]][[1]], legend = FALSE)

na_dat2 <- na_dat %>% drop_na(Phyto)
head(na_dat2) ##all variables give NA-values at same points (= check that cropping was done well)

rm(na_dat, na_dat2)

#remove observations without corresponding environmental values
full_dat_sp <- full_dat_sp %>%
  drop_na()

#3. Sampling bias - filtering -----
##3.1 adult data ----
### sampling bias: spatial filtering ----
species <- tibble(scientific = c("Clupea harengus", "Scomber scombrus", "Alosa fallax", "Dicentrarchus labrax"),
                  simple = c("herring", "mackerel", "twaite_shad", "seabass"))

#thin towards minimum distance of 10 NM or 18.52 km
presences_list <- list()
for (s in 1:nrow(species)) {
  tmp_df <- full_dat %>% filter(scientificname == species$scientific[s])
  presences_list[[s]] <- spThin::thin(tmp_df, lat.col = "lat", long.col = "lon", 
                                      spec.col = "scientificname", thin.par = 18.52, reps = 1,
                                      write.files = FALSE, locs.thinned.list.return = TRUE)[[1]]
  names(presences_list)[s] <- species$simple[s]
  print(s)
}

mapview(presences_list[[1]]$Longitude, presences_list[[1]]$Latitude, crs = 'epsg:4326')
mapview(presences_list[[2]]$Longitude, presences_list[[2]]$Latitude, crs = 'epsg:4326')
mapview(presences_list[[3]]$Longitude, presences_list[[3]]$Latitude, crs = 'epsg:4326')
mapview(presences_list[[4]]$Longitude, presences_list[[4]]$Latitude, crs = 'epsg:4326')

### sampling bias: environmental filtering ----
table(full_dat$scientificname)
# Alosa fallax      Clupea harengus Dicentrarchus labrax     Scomber scombrus 
#          546                12898                 1043                 7275 

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
    select(any_of(v_list)) %>%
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
# Alosa fallax      Clupea harengus Dicentrarchus labrax     Scomber scombrus 
#          160                  400                  160                  400

table(reduced_dat$scientificname, reduced_dat$Month)
#                        1   2   3   7   8   9  10  11  12
# Alosa fallax          32  71   0   4  24   7  15   7   0
# Clupea harengus       60  98  31  29  85  21  30  32  14
# Dicentrarchus labrax  15  17  25  10   5  16  50  17   5
# Scomber scombrus      23  46  18  40 122  37  31  61  22

##3.2 spawning data ----
### sampling bias: spatial filtering ----
species_sp <- tibble(scientific = c("Clupea harengus", "Scomber scombrus"),
                  simple = c("herring", "mackerel"))

#thin towards minimum distance of 10 NM or 18.52 km
presences_list_sp <- list()
for (s in 1:nrow(species_sp)) {
  tmp_df <- full_dat_sp %>% filter(scientificname == species_sp$scientific[s])
  presences_list_sp[[s]] <- spThin::thin(tmp_df, lat.col = "lat", long.col = "lon", 
                                      spec.col = "scientificname", thin.par = 18.52, reps = 1,
                                      write.files = FALSE, locs.thinned.list.return = TRUE)[[1]]
  names(presences_list_sp)[s] <- species_sp$simple[s]
  print(s)
}

mapview(presences_list_sp[[1]]$Longitude, presences_list_sp[[1]]$Latitude, crs = 'epsg:4326')
mapview(presences_list_sp[[2]]$Longitude, presences_list_sp[[2]]$Latitude, crs = 'epsg:4326')

### sampling bias: environmental filtering ----
table(full_dat_sp$scientificname)
# Clupea harengus Scomber scombrus 
#            5242             6020 

red_dat_sp <- list()
for (s in 1:nrow(species_sp)) {
  name <- species_sp$simple[s]
  
  tmp_full_dat <- left_join(presences_list_sp[[which(names(presences_list_sp) == species_sp$simple[s])]] %>% 
                              mutate(lon = Longitude, lat = Latitude), full_dat_sp, by = c("lon","lat")) %>%
    filter(scientificname == species_sp$scientific[s])
  
  df_pr <- tmp_full_dat[0,]
  df_rem <- tmp_full_dat
  
  #remove variables that have only one unique value
  tmp_df <- tmp_full_dat %>% 
    select(any_of(v_list)) %>%
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
  while (nrow(df_pr) < 400) {
    mal_dis <- mahalanobis.dist(df_pr[,12:16],df_rem[,12:16])
    mal_dis_tot <- colSums(mal_dis)
    id <- as.vector(which(mal_dis_tot == max(mal_dis_tot), useNames = FALSE))
    if (length(id) > 1) {id <- sample(id,1)}
    df_pr <- rbind(df_pr,df_rem[id,])
    df_rem <- df_rem[-id,]
  }
  red_dat_sp[[s]] <- df_pr
  print(name)
}

reduced_dat_sp <- bind_rows(red_dat_sp)
table(reduced_dat_sp$scientificname)
# Clupea harengus Scomber scombrus 
#             400              400

table(reduced_dat_sp$scientificname, reduced_dat_sp$Month)
#                    1   2   3   4   5   6   7   9  10  12
# Clupea harengus  145   0   0   0   0   0   0 216   3  36
# Scomber scombrus   0   4  46  43 120 158  29   0   0   0

#4. Save data ----
# stack list, mask and ICES areas shapefile
save(st_list_NEA_cl,
     final_common_mask, ices_shp,  
     file = "SAVE/final_st_list.Rdata")

# biological data before and after environmental filtering
save(presences_list_sp, presences_list, file = "SAVE/presences_after_spatial_filtering.Rdata")
save(full_dat, reduced_dat, file = "SAVE/bio_dat_adults.Rdata")
save(full_dat_sp, reduced_dat_sp, file = "SAVE/bio_dat_larvae.Rdata")
