pckgs <- c("raster","sp","proj4","ncdf4","car","mgcv","dismo","rJava","ENMeval","randomForest", "VSURF",
           "randomForestSRC","ggRandomForests","ggplot2","boot","gstat","mgcv","SDMtune","ENMeval",
           "tidyverse","mapdata","base","tidync", "sf", "usdm", "mapview", "tictoc", "ape", "spdep",
           "spThin", "StatMatch")
lapply(pckgs, require, character.only = TRUE)
rm(pckgs)

#1. Read data -----

## Biological data ----
#adults
setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Scripts data-driven approach/Herring/Scripts/Preprocessing/fishdish functions pimped/INPUT")
load("her_mack_tws_seab_adults_processed.Rdata")

#spawning
setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/GIS/Herring/biotic layers/OBIS on herring")
herring_sp <- read.csv("OBIS_larvae_herring_2000-2020.csv")

setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/GIS/Mackerel/biotic layers/OBIS on mackerel")
mackerel_sp <- read.csv("OBIS_larvae_mackerel_2000-2020.csv")

dati_sp <- rbind(herring_sp %>% dplyr::select(lon, lat, Year = year, Month = month, Day = day) %>%
                   mutate(scientificname = "Clupea harengus"),
                 mackerel_sp %>% dplyr::select(lon, lat, Year = year, Month = month, Day = day) %>%
                   mutate(scientificname = "Scomber scombrus"))
rm(herring_sp, mackerel_sp)

###plot data and remove suspicious points ----
#1. general
x <- data.frame(lon = dati_ad$lon, 
                lat = dati_ad$lat, 
                species = dati_ad$scientificname)
coordinates(x) <- ~lon+lat
crs(x) <- "epsg:4326"
mapview(x, zcol = "species")
#remove Baltic Sea observations (min lon 9.3656° according to marine regions)
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


#spawning
table(dati_sp$scientificname)
mapview(dati_sp %>% filter(scientificname == "Clupea harengus") %>% select(lon) %>% pull,
        dati_sp %>% filter(scientificname == "Clupea harengus") %>% select(lat) %>% pull,
        crs = "epsg:4326")
mapview(dati_sp %>% filter(scientificname == "Scomber scombrus") %>% select(lon) %>% pull,
        dati_sp %>% filter(scientificname == "Scomber scombrus") %>% select(lat) %>% pull,
        crs = "epsg:4326")


## Covariates ----

# only run once (creates stack list for BPNS and NEA)
# -- BPNS -- #
# results in: st_list_bpns[[m]][[v]] with m for month and v for variable
# setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Scripts/Herring/Input-output_files Herring/1.Input layers/")
# fls_bpns <- tibble(full_names    =   list.files(pattern = ".tif$|.grd$"),
#                    short_names   =   str_remove(full_names, ".tif$|.grd$"),
#                    var           =   str_remove_all(short_names, "_\\d{4}|_\\d{2}|_\\d"),
#                    mon           =   as.numeric(str_remove(str_extract(short_names, "_\\d{2}$|_\\d$"),"_")),
#                    year          =   as.numeric(str_remove(str_extract(short_names, "_\\d{4}"),"_")))
# 
# 
# st_list_bpns <- list()
# for (m in 1:12) {
#   tmp_fls <- fls_bpns %>%
#     filter(mon == m | is.na(mon))
#   tmp_st <- stack(tmp_fls$full_names)
#   #for now, use only most-likely scenario
#   names(tmp_st)[which(str_detect(names(tmp_st), "max_SSV"))] <- "mSSV"
#   tmp_st <- tmp_st[[-which(str_detect(names(tmp_st),"min|max"))]]
#   names(tmp_st) <- tmp_fls$var
#   st_list_bpns[[m]] <- tmp_st
# }
# rm(tmp_st, tmp_fls, fls_bpns)



# 
# -- North-East Atlantic -- #
# results in: st_list_NEA[[y]][[m]][[v]] with y for year, m for month and v for variable
# setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Scripts data-driven approach/Herring/Input-output_files/1.Resampled files")
# 
# fls_NEA <- tibble(full_name = list.files(pattern = ".grd$"),
#                   shrt_name = str_remove_all(full_name, ".grd$"))
# 
# st_list_NEA <- list()
# tmp_st_list <- list()
# for (y in 1:length(2000:2020)) {
#   for (m in 1:12) {
#     tmp_st_m <- stack()
#     for (v in 1:nrow(fls_NEA)) {
#       tmp_st <- stack(fls_NEA$full_name[v])
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
#   print(paste("processed", y))
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
# setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Scripts data-driven approach/Herring/Input-output_files/1.Resampled files")
# save(st_list_NEA, st_list_bpns, file = "stacks_covariates_not_cropped.rData")

# load files directly
setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Scripts data-driven approach/Herring/Input-output_files/1.Resampled files")
load("stacks_covariates_not_cropped.rData")

# #create average and sd layers for each environmental variable
# for (v in 1:nrow(fls_NEA)) {
#   tmp_st <- stack(fls_NEA$full_name[v])
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
#                   paste0("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Reporting/extra_GIS_layers/average_sd_vars/NEA_average_",
#                          fls_NEA$shrt_name[v],"_",
#                          m,".tif"), overwrite = T)
#       writeRaster(sd_NEA, 
#                   paste0("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Reporting/extra_GIS_layers/average_sd_vars/NEA_sd_",
#                          fls_NEA$shrt_name[v],"_",
#                          m,".tif"), overwrite = T)
#     }
#   }
# }


### crop layers ----

#### BPNS ----
plot(st_list_bpns[[1]]) #needs cropping
bpns_shp <- st_read("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/GIS/Herring/other/eez/eez.shp")
st_list_bpns2 <- lapply(st_list_bpns, mask, mask = bpns_shp)
plot(st_list_bpns2[[1]])


#### NEA ----
plot(st_list_NEA[[1]][[1]])

ices_shp <- st_read("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/GIS/Herring/other/ICES areas/ICES_Areas_20160601_cut_dense_3857.shp")
mapview(ices_shp)

ices_selection_ad <- ices_shp[which(ices_shp$Area_27 %in% dati_ad$Area_27),]
ices_selection_ad <- ices_selection_ad %>%
  st_transform("epsg:4326") %>%
  st_sf %>%
  st_cast
mapview(ices_selection_ad)


#add ICES area 27 to spawning data
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

## crop raster towards the ICES regions where observations are present
## also crop to common non-NA values of all rasters (common mask)
## TODO: move to resampling script!
## change to also add spawning ones?

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
# check categorical variables seperately (different test)
# long loop, only run once
df_ene <- as.data.frame(st_list_NEA2[[1]][[1]][[8]]) %>%
  rename(labels = seabed_energy_2000_1) %>%
  left_join(seabed_ene, by = "labels") %>%
  select(-labels) %>%
  dplyr::pull()
df_sub <- as.data.frame(st_list_NEA2[[1]][[1]][[9]]) %>%
  rename(labels = seabed_substrate_2000_1) %>%
  left_join(seabed_hab, by = "labels") %>%
  select(-labels) %>%
  dplyr::pull()

for (y in 1:length(2000:2020)) {
  for (m in 1:12) {
    for (v in (1:length(names(st_list_NEA2[[1]][[1]])))[-c(8,9)]){
      #skip loop if layer has only 1 unique value
      if(length(unique(na.omit(pull(as.data.frame(st_list_NEA2[[y]][[m]][[v]]))))) == 1) next
      
      #Kruskal Wallis test to test for correlation
      t_s <- kruskal.test(df_sub, pull(as.data.frame(st_list_NEA2[[y]][[m]][[v]])))
      t_e <- kruskal.test(df_ene, pull(as.data.frame(st_list_NEA2[[y]][[m]][[v]])))
      
      if(t_s$p.value < 0.05) {
        print(paste("substrate ~", names(st_list_NEA2[[1]][[1]])[v], (2000:2020)[y], m))
        print(t_s)
      }
      if(t_e$p.value < 0.05) {
        print(paste("energy ~", names(st_list_NEA2[[1]][[1]])[v], (2000:2020)[y], m))
        print(t_e)
      }
    }
  }
  print(y)
}

#don't include seabed energy & phyto/chlorophyll in the same model



#### per month ----
#function
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

# chlorophyll and phytoplankton are highly correlated (as expected) --> remove chlorophyll, 
# since we have response curve for phytoplankton in mechanistic approach
retained_vars <- all_vars[-which(all_vars == "Phyto")]
vif_t2 <- f_vif(st_list_NEA2, retained_vars)
any(vif_t2[,-c(1,2)] > 10)
vif_t2$Var[apply(vif_t2[,-c(1,2)] > 10, 1, any)]
vif_t2[,1:3]

# temperature and dissolved oxygen are highly correlated --> remove dissolved oxygen, 
# since temperature is an important variable to keep
retained_vars <- all_vars[-which(all_vars == "Phyto" | all_vars == "DO" )]
vif_t3 <- f_vif(st_list_NEA2, retained_vars)
vif_t3[,1:3]

# remove Euphotic depth
retained_vars <- all_vars[-which(all_vars == "Phyto" | all_vars == "DO" | all_vars == "EuphD")]
vif_t3 <- f_vif(st_list_NEA2, retained_vars)
vif_t3[,1:3]

#conclusion: remove Phyto, DO and EuphD


#### using correlation plots ----
# Chl and depth

tib_all_v <- tibble(a = rep(1,3846528)) %>% select()
for (v in 1:length(all_vars)) {
  nm <- all_vars[v]
  tmp_x <- NULL
  for (y in 1:length(2000:2020)) {
    for (m in 1:12) {
      tmp_st <-  st_list_NEA2[[y]][[m]]
      tmp_x <- append(tmp_x, values(tmp_st[[which(str_detect(names(tmp_st), nm))]]))
    }
  }
  tib_all_v <- tib_all_v %>% add_column(tmp_x)
  colnames(tib_all_v)[v] <- all_vars[v]
  print(nm)
}

tib_all_v$Month <- rep(1:12, times = 21, each = ncell(st_list_NEA2[[1]][[1]])) 
tib_all_v$Year <- rep(2000:2020, times = 1, each = ncell(st_list_NEA2[[1]][[1]])*12) 

tib_all_v <- tib_all_v %>% na.omit %>% as.data.frame()
head(tib_all_v)


vif(tib_all_v)
s <- sample(nrow(tib_all_v), 10000)
pairs(tib_all_v %>% select(Chl, Phyto) %>% filter(row_number() %in% s))

#remove SSS
vif(tib_all_v %>% select(-Phyto))
cor(tib_all_v %>% select(-Phyto))

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(tib_all_v %>% filter(row_number() %in% s) %>% select(-Chl), 
      lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE)
cor(tib_all_v)
vifcor(tib_all_v)

#conclusion: remove Phyto


# # say we remove distance from coast
# retained_vars <- all_vars[-which(all_vars == "Phyto" | all_vars == "DO" | all_vars == "distance_from_coast" )]
# vif_t5 <- f_vif(st_list_NEA2, retained_vars)
# any(vif_t5[,-c(1,2)] > 10)
# tibble(var = retained_vars,
#        max = apply(vif_t5[,-c(1,2)], 1, max, na.rm = T))

# All_50_SPDF <- as(stack_50, "SpatialPixelsDataFrame")
# All_50_DF <- as.data.frame(All_50_SPDF)
# head(All_50_DF)

###conclusion: remove Phyto from stack list ----

v_list <- all_vars[-which(all_vars == "Phyto")]

st_list_NEA_cl <- st_list_NEA2
for (y in 1:length(2000:2020)) {
  for (m in 1:12) {
    st_list_NEA_cl[[y]][[m]] <- st_list_NEA2[[y]][[m]][[which(str_detect(names(st_list_NEA2[[y]][[m]]), paste(v_list, collapse = "|")))]]
  }
}

#2. sample covariates with biological data ----

##2.1 adult data ----
## remove duplicates ----
dati_ad <- dati_ad[-which(duplicated(dati_ad %>% select(Year, Month, lon, lat, scientificname))),]

## sampling ----
coords_dat <- dati_ad %>% select(lon, lat, Year, Month, Day, scientificname, Gear, Ship, Area_27)
full_dat <- data.frame()
for (y in 1:length(2000:2020)) {
  for (m in 1:12) {
    tmp_coords_dat <- coords_dat %>%
      filter(Year == c(2000:2020)[y],
             Month == m)
    tmp_r <- st_list_NEA2[[y]][[m]]
    tmp_extract <- raster::extract(st_list_NEA2[[y]][[m]], 
                                   as.data.frame(tmp_coords_dat[,c(1,2)]), method = "simple")
    colnames(tmp_extract) <- colnames(tmp_extract) %>%
      str_remove_all("_\\d{4}_\\d{2}|_\\d{4}_\\d")
    tmp_df <- data.frame(tmp_coords_dat, tmp_extract)
    full_dat <- rbind(full_dat, tmp_df)
  }
  print(y)
}
rm(tmp_extract, tmp_coords_dat, tmp_df)

## NA values ----
any(is.na(full_dat$Gear)) #no NA values
any(is.na(full_dat$Ship)) #no NA values
na_dat <- full_dat %>% filter(is.na(depth) | is.na(distance_from_coast) | 
                                is.na(EuphD) | is.na(max_SSV) | is.na(Chl) | 
                                is.na(SSS) | is.na(SST) | is.na(ZooPl))

head(na_dat)
mapview::mapview(na_dat$lon, na_dat$lat, crs = CRS("epsg:4326"), legend = FALSE) +
  mapview::mapview(st_list_NEA_cl[[1]][[3]][[1]], legend = FALSE)
mapview::mapview(na_dat$lon, na_dat$lat, crs = CRS("epsg:4326"), legend = FALSE) +
  mapview::mapview(st_list_NEA_cl[[1]][[1]][[3]], legend = FALSE)

na_dat2 <- na_dat %>% drop_na(Phyto)
head(na_dat2) ##all variables give na-values at same points (check that cropping was done well)

rm(na_dat, na_dat2)

full_dat <- full_dat %>%
  drop_na()

##2.2 spawning data ----
## remove duplicates ----
dati_sp <- dati_sp[-which(duplicated(dati_sp %>% select(Year, Month, lon, lat, scientificname))),]

## sampling ----
coords_dat_sp <- dati_sp %>% select(lon, lat, Year, Month, Day, scientificname, Area_27)
full_dat_sp <- data.frame()
for (y in 1:length(2000:2020)) {
  for (m in 1:12) {
    tmp_coords_dat <- coords_dat_sp %>%
      filter(Year == c(2000:2020)[y],
             Month == m)
    tmp_r <- st_list_NEA2[[y]][[m]]
    tmp_extract <- raster::extract(st_list_NEA2[[y]][[m]], 
                                   as.data.frame(tmp_coords_dat[,c(1,2)]), method = "simple")
    colnames(tmp_extract) <- colnames(tmp_extract) %>%
      str_remove_all("_\\d{4}_\\d{2}|_\\d{4}_\\d")
    tmp_df <- data.frame(tmp_coords_dat, tmp_extract)
    full_dat_sp <- rbind(full_dat_sp, tmp_df)
  }
  print(y)
}
rm(tmp_extract, tmp_coords_dat, tmp_df)

## NA values ----
na_dat <- full_dat_sp %>% filter(is.na(depth) | is.na(distance_from_coast) | is.na(DO) | 
                                is.na(EuphD) | is.na(max_SSV) | is.na(Chl) | is.na(Phyto) | 
                                is.na(seabed_energy) | is.na(seabed_substrate) | 
                                is.na(SSS) | is.na(SST) | is.na(ZooPl))

head(na_dat)
mapview::mapview(na_dat$lon, na_dat$lat, crs = CRS("epsg:4326"), legend = FALSE) +
  mapview::mapview(st_list_NEA_cl[[1]][[3]][[1]], legend = FALSE)
mapview::mapview(na_dat$lon, na_dat$lat, crs = CRS("epsg:4326"), legend = FALSE) +
  mapview::mapview(st_list_NEA_cl[[1]][[1]][[3]], legend = FALSE)

na_dat2 <- na_dat %>% drop_na(Phyto)
head(na_dat2) ##all variables give na-values at same points (check that cropping was done well)

rm(na_dat, na_dat2)

full_dat_sp <- full_dat_sp %>%
  drop_na()

#3. Sampling bias - filtering -----
##3.1 adult data ----
#first remove months that are not represented in surveys (quarter 2, months 4-6)
table(full_dat$Month)
#    1    2    3    4    6    7    8    9   10   11   12 
# 1937 4317 1887   64    5 1314 5468 1549 2282 2512  508 

full_dat <- full_dat %>% filter(!(Month %in% c(4,6)))

table(full_dat$Month)
#    1    2    3    7    8    9   10   11   12 
# 1937 4317 1887 1314 5468 1549 2282 2512  508 

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


## outliers in species observations ----
# do not run

# f_IQ <- function(var) {
#   qs <- quantile(var)
#   IQD <- qs[4] - qs[2]   #interquantile distance equals 75% minus 25% interval
#   thresholds <- c(qs[2] - IQD*3, qs[4] + IQD*3)
#   c(which(var < thresholds[1]), which(var > thresholds[2]))
# }
# 
# plot_list <- list()
# for (s in 1:nrow(species)) {
#   tmp_full_dat <- left_join(presences_list[[which(names(presences_list) == species$simple[s])]] %>% 
#                               mutate(lon = Longitude, lat = Latitude), full_dat, by = c("lon","lat")) %>%
#     filter(scientificname == species$scientific[s])
#   
#   reduced_dat <- red_dat[[s]][0,]
#   red_dat_sp <- tmp_full_dat[0,]
#   tmp_full_dat <- red_dat[[s]]
# 
#   num_dat <- tmp_full_dat %>% select(where(is.numeric)) %>%  #select all numeric variables, but not lon, lat, longitude and latitude and month
#     select(-lon, -lat, -Longitude, -Latitude, -Month, -windfarms, -seabed_energy, -seabed_substrate)
#   
#   rm_d <- unlist(apply(num_dat, 2, f_IQ))
#   rm_d_tib <- tibble(variable = str_remove(names(rm_d), ".\\d*$"),
#                      rnum = unname(rm_d),
#                      val = rep(1))
# 
#   plot_dat <- num_dat %>% mutate(rnum = row_number()) %>% 
#     gather(key = "variable",
#            value = "value",
#            -rnum) %>%
#     left_join(rm_d_tib, by = c("rnum","variable")) %>%
#     mutate(val = ifelse(is.na(val), "in", "out"))
# 
#   plot_list[[s]] <- ggplot(plot_dat, aes(x = value, y = rnum, colour = val)) +
#     geom_point(size = 1) +
#     theme_bw() +
#     theme(
#       panel.grid.major.x = element_blank(),
#       panel.grid.minor.x = element_blank() 
#     ) +
#     facet_wrap(~variable, scales = "free") +
#     labs(y= "Order of the data", x = "Value of the variable")
# 
#   tmp_full_dat <- tmp_full_dat[-rm_d_tib$rnum,]
#   reduced_dat <- rbind(reduced_dat, tmp_full_dat)
# }
# 
# table(reduced_dat$scientificname)
# table(reduced_dat$scientificname, reduced_dat$Month)
# reduced_dat <- reduced_dat %>% filter(!Month %in% c(4,6))
# nrow(reduced_dat)
# reduced_dat <- reduced_dat %>% 
#   filter(!(Month %in% c(3,7,9,11,12) & scientificname == "Alosa fallax"),
#          !(Month %in% c(7,8,12) & scientificname == "Dicentrarchus labrax"))
# nrow(reduced_dat)
# table(reduced_dat$scientificname, reduced_dat$Month)
# 
# windows()
# plot_list[[1]]
# plot_list[[2]]
# plot_list[[3]]
# plot_list[[4]]
# 
# c(sum(nrow(presences_list[[1]]), nrow(presences_list[[2]]), 
#       nrow(presences_list[[3]]), nrow(presences_list[[4]])), 
#   nrow(reduced_dat))
# #removed 100 outliers in total

### sampling bias: environmental filtering ----
table(full_dat$scientificname)
# Alosa fallax      Clupea harengus Dicentrarchus labrax     Scomber scombrus 
# 558               12949           1050                     7286 

num_pr_vector <- c(350,350,160,160)

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
#          160                  200                  160                  200

table(reduced_dat$scientificname, reduced_dat$Month)
#  1   2   3   7   8   9  10  11  12 
# 84 175  58  43 107  56 104  72  21  

setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Scripts data-driven approach/Herring/Input-output_files/1.Resampled files")

## save data ----
save(st_list_NEA2,  
     file = "final_st_lists.Rdata")

save(st_list_bpns2, full_dat, reduced_dat, ices_selection_ad,
     final_common_mask, ices_shp, bpns_shp, 
     file = "bio_dat_co.Rdata")


# OUTLIERS IN X Y

# difficult to check because of huge size of rasters --> TBS
# also look at plots from Zuur!
# par(mfrow = c(3,3))
# 
# boxplot(st_list_NEA_cl[[1]][[1]]$depth_2000_1,ylab='depth')
# boxplot(stack_500$slope,data=stack_500,ylab='slope')
# boxplot(stack_500$bpi_fine,data=stack_500,ylab='bpi_fine')
# boxplot(stack_500$bpi_broad,data=stack_500,ylab='bpi_broad')
# boxplot(stack_500$eastness,data=stack_500,ylab='eastness')
# boxplot(stack_500$northness,data=stack_500,ylab='northness')
# boxplot(stack_500$roughness,data=stack_500,ylab='roughness')
# 
# windows()
# dotchart(as.matrix(stack_500$slope), xlab = "slope",ylab = "Order of the data")
# dotchart(as.matrix(stack_500$bpi_fine), xlab = "bpi_fine",ylab = "Order of the data")
# dotchart(as.matrix(stack_500$roughness), xlab = "roughness",ylab = "Order of the data")
# 
# # homogeniety in Y
# par(mfrow=c(3,3))
# plot(isi2_500$depth,isi2_500$PresAbs, ylab = "PresAbs", xlab = "depth")
# plot(isi2_500$slope,isi2_500$PresAbs, ylab = "PresAbs", xlab = "slope")
# plot(isi2_500$bpi_fine,isi2_500$PresAbs, ylab = "PresAbs", xlab = "bpi_fine")
# plot(isi2_500$bpi_broad,isi2_500$PresAbs, ylab = "PresAbs", xlab = "bpi_broad")
# plot(isi2_500$eastness,isi2_500$PresAbs, ylab = "PresAbs", xlab = "eastness")
# plot(isi2_500$northness,isi2_500$PresAbs, ylab = "PresAbs", xlab = "northness")
# plot(isi2_500$roughness,isi2_500$PresAbs, ylab = "PresAbs", xlab = "roughness")
# plot(isi2_500$vms2009_2014,isi2_500$PresAbs, ylab = "PresAbs", xlab = "vms2009_2014")
# plot(isi2_500$vms2015_2018,isi2_500$PresAbs, ylab = "PresAbs", xlab = "vms2015_2018")
# 
# plot(isi2_500$vms2015_2018,isi2_500$NumberColonies, ylab = "Number of colonies", xlab = "vms2015_2018")
# plot(isi2_500$vms2015_2018,isi2_500$NumberColonies, ylab = "Number of colonies", xlab = "vms2015_2018")
# 
# par(mfrow=c(1,2))
# hist(isidella$NumberColonies)
# qqnorm(isidella$NumberColonies); qqline(isidella$NumberColonies)
# plot(table(isidella$NumberColonies), type = "h",xlab = "Observed values", ylab = "Frequency")
# 
# #spatial independence of Y
# par(mfrow=c(1,3))
# plot(NumberColonies~ Lon_UTM + Lat_UTM, data=isi2_500, type='l',lwd=2)
# library(gstat)
# 
# spat.isi2_500 <- data.frame(isi2_500$PresAbs , isi2_500$Lon_UTM, isi2_500$Lat_UTM)
# coordinates(spat.isi2_500)<-c('isi2_500.Lon_UTM','isi2_500.Lat_UTM')
# 
# par(mfrow = c(1,2))
# vario1<-variogram(isi2_500.PresAbs~1,data=spat.isi2_500)
# plot(vario1)
# vario2<-variogram(isi2_500.PresAbs~1,data=spat.isi2_500, alpha=c(0,45,90,135))
# plot(vario2)

#create sample for EDITO guys
# write.csv((full_dat %>% select(-scientificname))[c(1:200),], "C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Presentations/Other/20230917-meeting with IT/data_extract.csv")



##3.2 spawning data ----
#first remove months that are not represented in surveys (month 8)
table(full_dat_sp$scientificname, full_dat_sp$Month)
#                      1     2     3     4     5     6     7     8     9    10    12
# Clupea harengus   2253     0     0     0     0     0     0     0  4420    95  1012
# Scomber scombrus     0   299  3255  7550 10840 10852  1770     1     0     0     0 

full_dat_sp <- full_dat_sp %>% filter(!(Month %in% c(8)))

table(full_dat_sp$scientificname, full_dat_sp$Month)
#                      1     2     3     4     5     6     7     9    10    12
# Clupea harengus   2253     0     0     0     0     0     0  4420    95  1012
# Scomber scombrus     0   299  3255  7550 10840 10852  1770     0     0     0


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

setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Scripts data-driven approach/Herring/Input-output_files/1.Resampled files")
save(presences_list_sp, presences_list, file = "presences_after_spatial_filtering.rData")

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
  while (nrow(df_pr) < 350) {
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
#             350              350

table(reduced_dat_sp$scientificname, reduced_dat_sp$Month)
#                    1   2   3   4   5   6   7   9  10  12
# Clupea harengus  123   0   0   0   0   0   0 171   7  49
# Scomber scombrus   0   4  41  42 106 150   7   0   0   0

setwd("C:/Users/ward.standaert/OneDrive - VLIZ/BAR ecological modelling/Scripts data-driven approach/Herring/Input-output_files/1.Resampled files")

## save data ----
save(full_dat_sp, reduced_dat_sp, file = "bio_dat_co_spawning.Rdata")
