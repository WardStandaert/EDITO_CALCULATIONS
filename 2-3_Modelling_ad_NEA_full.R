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

#0. load data & functions ----

load("~/INPUT/final_st_lists.Rdata")
load("~/INPUT/bio_dat_co_sp&env_sampling.Rdata")

setwd("~/OUTPUT/all_data_model/adult")

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
species <- tibble(scientific = c("Clupea harengus", "Scomber scombrus", "Alosa fallax", "Dicentrarchus labrax"),
                  simple = c("herring", "mackerel", "twaite_shad", "seabass"))

# Variables selected per species
# Twaite shad includes all vars except for phytoplankton because the variable is correlated with seabed_energy

v_list <- list()
v_list$herring <- c("SST", "SSS", "windfarms", "ZooPl", "Phyto", "seabed_energy", "seabed_substrate", "depth")
v_list$mackerel <- c("SST", "SSS", "windfarms", "max_SSV", "EuphD", "ZooPl", "depth")
v_list$twaite_shad <- c("Quarter_sin", "Quarter_cos", "SST", "SSS", "windfarms", "ZooPl", "seabed_energy", "seabed_substrate", "depth") 
v_list$seabass <- c("SST", "SSS", "windfarms", "ZooPl", "depth")


# Circular encoding for month where january (1) and december (12) are similar
full_dat <- full_dat %>%
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

#1. extract environmental values for presences ---- 

#convert month to factor

pr_pa <- list()
pr_coord <- list()
pr_predv <- list()
pr_pa_coord_predv <- list()

for (s in 1:nrow(species)) {
  tmp_pres  <- full_dat %>%
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
  tmp_dat <- full_dat %>%
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

# start from SAVE
# save(pr_pa, pr_coord, pr_predv,
#      bg_pa, bg_coord, bg_predv,
#      file = "SAVE/NEA_presence_background.Rdata")

# SAVE ----
load("SAVE/NEA_presence_background.Rdata")

#3 ENMeval model creation  ----
#convert month to character bc algorithm cant handle factors that look like numeric
eval_res_list <- list()
tmp_pr_predv <- list()
tmp_bg_predv <- list()
for (s in 1) {
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

print(eval_res_list[[1]] %>%  eval.results() %>% filter(delta.AICc == 0))

#5. Spatial autocorrelation ----

## presence data ----
morans_I_list <- list()
for (s in 1) {
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


morans_I_list[[1]]   
#Before filtering:          Moran's I value of 0.18, expected -0.00008; significant
# morans_I_list[[2]]   
# #Before filtering:          Moran's I value of 0.16, expected -0.0001; significant
# #After sp & geo filtering:  Moran's I value of 0.20, expected -0.01; significant
# morans_I_list[[3]]   
# #Before filtering:          Moran's I value of 0.34, expected -0.002; significant
# #todo After sp & geo filtering:  Moran's I value of 0.20, expected -0.01; significant
# morans_I_list[[4]]   
# #Before filtering:          Moran's I value of 0.48, expected -0.001; significant
# #After sp & geo filtering:  Moran's I value of 0.22, expected -0.01; significant


