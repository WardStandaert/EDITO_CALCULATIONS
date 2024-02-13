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
setwd("DATA/")
load("SAVE/final_st_list.Rdata")
load("SAVE/bio_dat_adults.Rdata")
load("SAVE/presences_after_spatial_filtering.rData")
rm(presences_list_sp)

# Variable names
v_list_all <- str_remove_all(names(st_list_NEA_cl[[1]][[1]]) , "_\\d{4}_\\d{2}|_\\d{4}_\\d")
# Species names
species <- tibble(scientific = c("Clupea harengus", "Scomber scombrus" , "Alosa fallax", "Dicentrarchus labrax"),
                  simple = c("herring", "mackerel" , "twaite_shad", "seabass"))

# Variables selected per species through literature review
v_list <- list()
v_list$herring <- c("SST", "SSS", "windfarms", "ZooPl", "Phyto", "seabed_energy", "seabed_substrate", "depth")
v_list$mackerel <- c("SST", "SSS", "windfarms", "max_SSV", "ZooPl", "depth")
v_list$twaite_shad <- c("SST", "SSS", "windfarms", "ZooPl", "seabed_energy", "seabed_substrate", "depth")
v_list$seabass <- c("SST", "SSS", "windfarms", "ZooPl", "depth")

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

pr_pa <- list()
pr_coord <- list()
pr_predv <- list()

for (s in 1:nrow(species)) {
  tmp_pres  <- reduced_dat %>%
    filter(scientificname == species$scientific[s])
  
  pr_pa[[s]] <- rep(1, nrow(tmp_pres))
  pr_coord[[s]] <- tmp_pres %>% dplyr::select(lon, lat, Year, Month)
  pr_predv[[s]] <- tmp_pres %>% dplyr::select(all_of(v_list[[which(names(v_list) == species$simple[s])]]))
  
  names(pr_pa)[s] <- names(pr_coord)[s] <- names(pr_predv)[s] <- species$simple[s]
}

# check the outcome
lapply(pr_pa, head, n=2)              #presence points presence/absence values per species
lapply(pr_coord, head, n=2)           #presence points coordinates per species
lapply(pr_predv, head, n=2)           #presence points environmental variable extracts per species


#2. create background points and extract environmental values ----
set.seed(123)

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
                               Month = sample(unique(tmp_dat$Month), 10*length(pr_pa[[s]]), replace = TRUE))

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

# check the outcome
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

# start from SAVE
save(pr_pa, pr_coord, pr_predv,
     bg_pa, bg_coord, bg_predv,
     file = "SAVE/NEA_presence_background_adults.Rdata")

# SAVE ----
load("SAVE/NEA_presence_background_adults.Rdata")

#3 ENMeval model creation  ----
#convert month to character bc algorithm cant handle factors that look like numeric
eval_res_list <- list()
tmp_pr_predv <- list()
tmp_bg_predv <- list()
for (s in 1:nrow(species)) {
  if(all(c("seabed_energy","seabed_substrate") %in% names(pr_predv[[s]]))) {
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
  
  #test different model settings using ENMeval
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

save(eval_res_list, file = "SAVE/NEA_ENMeval_outcomes_adults.Rdata")

# SAVE ----
load("SAVE/NEA_ENMeval_outcomes_adults.Rdata")

for(s in 1:nrow(species)) print(eval_res_list[[s]] %>%  eval.results() %>% filter(delta.AICc == 0))


#4 model evaluation AUC & TSS ----
set.seed(123)

AUC_maxent <- list()
TSS_maxent <- list()

set.seed(123)

for (s in 1:nrow(species)) {
  tic(species$simple[s])
  Presences <- pr_predv[[s]]
  Background <- bg_predv[[s]]
  
  if(all(c("seabed_energy","seabed_substrate") %in% names(pr_predv[[s]]))) {
    Presences <- pr_predv[[s]] %>%
      left_join(energy_lvl, by = "seabed_energy") %>%
      left_join(substr_lvl, by = "seabed_substrate") %>%
      dplyr::select(-seabed_energy, -seabed_substrate)
    
    Background <- bg_predv[[s]] %>%
      left_join(energy_lvl, by = "seabed_energy") %>%
      left_join(substr_lvl, by = "seabed_substrate") %>%
      dplyr::select(-seabed_energy, -seabed_substrate)
  }
  else {
    Presences <- pr_predv[[s]]
    Background <- bg_predv[[s]]
  }

  tmp_AUC <- list()
  tmp_TSS <- list()
  for (i in 1:10) {
    tic(paste("iteration", i))
    #Generate "k" groups:
    k<-4 #division of the data will be 75% Vs 25%

    groups_pres<-kfold(Presences,k) #Kfold divide the data, assigning every row to one of the K groups randomly. 
    groups_abs<-kfold(Background,k)
    
    #Four groups will be used to generate the model and the rest of the point (one group) will be used to evaluate it:
    EvalBg<-Background[groups_pres==1,]
    TrainBg<-Background[groups_pres!=1,]
    
    EvalPres<-Presences[groups_pres==1,]
    TrainPres<-Presences[groups_pres!=1,]

    #get model settings
    eval_res <- eval_res_list[[s]]
    opt.aicc <- eval.results(eval_res_list[[s]]) %>% filter(delta.AICc == 0)
    mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]    
    
    EvalBgRes  <- predict(mod, EvalBg, type = "cloglog")
    EvalPresRes  <- predict(mod, EvalPres, type = "cloglog")
    
    tmp_AUC[[i]]<-evaluate(c(EvalPresRes), c(EvalBgRes))
    umbral_maxent<-threshold( tmp_AUC[[i]])   
    tmp_TSS[[i]]<-evaluate(p=c(EvalPresRes), a=c(EvalBgRes), tr=umbral_maxent$spec_sens) 
    
    toc()
  }
  
  AUC_maxent[[s]] <- tmp_AUC
  TSS_maxent[[s]] <- tmp_TSS
  names(AUC_maxent)[s] <- names(TSS_maxent)[s] <- species$simple[s]
  toc()
}

AUC_maxent_vect   <- data.frame(matrix(nrow = 10, ncol = 4))
TPR_maxent_vect   <- data.frame(matrix(nrow = 10, ncol = 4))
TNR_maxent_vect   <- data.frame(matrix(nrow = 10, ncol = 4))
TSS_maxent_vect   <- data.frame(matrix(nrow = 10, ncol = 4))

for (s in 1:nrow(species)) {
  AUC_maxent_vect[,s]  <- sapply(AUC_maxent[[s]],function(x){slot(x,'auc')})
  TPR_maxent_vect[,s]  <- sapply(TSS_maxent[[s]],function(x){slot(x,'TPR')})
  TNR_maxent_vect[,s]  <- sapply(TSS_maxent[[s]],function(x){slot(x,'TNR')})
  TSS_maxent_vect[,s]  <- TPR_maxent_vect[,s] + TNR_maxent_vect[,s] - 1
  
  colnames(AUC_maxent_vect)[s] <- colnames(TPR_maxent_vect)[s] <- 
    colnames(TNR_maxent_vect)[s] <-  colnames(TSS_maxent_vect)[s] <- 
    species$simple[s]
}

sum_auc <- AUC_maxent_vect %>%
  gather(key = "variable",
         value = "value") %>%
  dplyr::group_by(variable) %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value))

sum_tss <- TSS_maxent_vect %>%
  gather(key = "variable",
         value = "value") %>% 
  dplyr::group_by(variable) %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value))

sum_TPR <- TPR_maxent_vect %>%
  gather(key = "variable",
         value = "value") %>% 
  dplyr::group_by(variable) %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value))

sum_TNR <- TNR_maxent_vect %>%
  gather(key = "variable",
         value = "value") %>% 
  dplyr::group_by(variable) %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value))

sum_all <- cbind(sum_auc, sum_tss[,c(2:3)], sum_TPR[,c(2:3)], sum_TNR[,c(2:3)])
colnames(sum_all) <- c("species",
                       "av_AUC", "sd_AUC", 
                       "av_TSS", "sd_TSS",
                       "av_TPR", "sd_TPR", 
                       "av_TNR", "sd_TNR")
sum_all[match(sum_all$species, species$simple),]

#save results
write.csv(sum_all, "3.MODEL_OUTPUT/ADULTS/validation_metrics.csv")


#5. Spatial autocorrelation ----
morans_I_list <- list()
for (s in 1:nrow(species)) {
  tic(paste0(species$simple[s], " - done"))
  eval_res <- eval_res_list[[s]]
  res <- eval.results(eval_res)
  opt.aicc <- res %>% filter(delta.AICc == 0)
  mod.seq <- eval.models(eval_res)[[opt.aicc$tune.args]]
  
  #2. make original prediction
  df_pred <- pr_predv[[s]]
  if(all(c("seabed_energy", "seabed_substrate") %in% names(df_pred))) {
      df_pred <- df_pred %>%
        mutate(sub_char = as.numeric(seabed_substrate)) %>%
        mutate(ene_char = as.numeric(seabed_energy)) %>%
        dplyr::select(-seabed_substrate, -seabed_energy)
      }
  
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
  
  toc()
}

morans_I_list[[1]]   
#Before filtering:          Moran's I value of 0.18, expected -0.00008; significant
#After sp & geo filtering:  Moran's I value of 0.09, expected -0.01; significant
morans_I_list[[2]]   
#Before filtering:          Moran's I value of 0.16, expected -0.0001; significant
#After sp & geo filtering:  Moran's I value of 0.22, expected -0.01; significant
morans_I_list[[3]]  
#Before filtering:          Moran's I value of 0.20, expected -0.002; significant
#After sp & geo filtering:  Moran's I value of 0.26, expected -0.01; significant
morans_I_list[[4]]  
#Before filtering:          Moran's I value of 0.21, expected -0.001; significant
#After sp & geo filtering:  Moran's I value of 0.21, expected -0.01; significant

#6. variable importance ----

#make big tibble with all possible values per variable

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
# tib_all_v$Year <- rep(2000:2020, times = 1, each = ncell(st_list_NEA_cl[[1]][[1]])*12)
# colnames(tib_all_v)[which(colnames(tib_all_v) %in% c("ene_char", "sub_char"))] <- c("seabed_energy", "seabed_substrate")
# 
# tib_all_v <- tib_all_v %>%
#   left_join(substr_lvl %>% mutate(seabed_substrate = as.numeric(seabed_substrate)), by = "seabed_substrate") %>%
#   left_join(energy_lvl %>% mutate(seabed_energy = as.numeric(seabed_energy)), by = "seabed_energy") %>%
#   dplyr::select(-seabed_substrate, -seabed_energy)
# 
# save(tib_all_v, file = "SAVE/tib_all_v.Rdata")
load("SAVE/tib_all_v.rData")

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
  if(all(c("seabed_energy", "seabed_substrate") %in% names(df_pred))) {
    df_pred <- df_pred %>%
      left_join(energy_lvl, by = "seabed_energy") %>%
      left_join(substr_lvl, by = "seabed_substrate") %>%
      dplyr::select(-seabed_substrate, -seabed_energy)
  }
  
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

# save(cor_maxnet, file = "SAVE/NEA_variable_importance_adults.Rdata")

# SAVE ----
load("SAVE/NEA_variable_importance_adults.Rdata")

#convert to percentages
var_imp_ad <- cor_maxnet[[1]]
var_imp_ad <- round(var_imp_ad/sum(var_imp_ad),2)
var_imp_ad
sum(var_imp_ad)

write.csv(var_imp_ad, "3.MODEL_OUTPUT/ADULTS/variable_importance.csv")

#7. Spatial predictions ----
for (m in 1:12) {
  pr_geo_y <- stack()
  for (y in 1:length(2000:2020)) {
    #get model
    eval_res <- eval_res_list[[1]]
    opt.aicc <- eval.results(eval_res_list[[1]]) %>% filter(delta.AICc == 0)
    mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]
    plot_st <- st_list_NEA_cl[[y]][[m]]
    
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
  
  # crop to extent of study area
  pr_geo_y <- crop(pr_geo_y, extent(-12, 10, 48, 62))
  
  # average & coefficient of variation NEA
  av_NEA <- stackApply(pr_geo_y, indices =  rep(1,nlayers(pr_geo_y)), mean, na.rm = T)
  sd_NEA <- stackApply(pr_geo_y, indices =  rep(1,nlayers(pr_geo_y)), sd, na.rm = T)
  uncert_NEA <- sd_NEA/av_NEA*100
  
  writeRaster(av_NEA, paste0("3.MODEL_OUTPUT/ADULTS/RASTERS/adult_herring_average_HSI_",m,".tif"), overwrite = T)
  writeRaster(sd_NEA, paste0("3.MODEL_OUTPUT/ADULTS/RASTERS/adult_herring_sd_HSI_",m,".tif"), overwrite = T)

  print(m)
}


#8. HSI variation BPNS ----
bpns_shp <- st_read("1.DOWNLOAD/shapefiles/BPNS/eez.shp", quiet = TRUE)
df_out <- data.frame(month = rep(month.abb, each = 9))

for (y in 1:length(2000:2020)) {
  pr_geo_y <- stack()
  for (m in 1:12) {
    #get model
    eval_res <- eval_res_list[[1]]
    opt.aicc <- eval.results(eval_res_list[[1]]) %>% filter(delta.AICc == 0)
    mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]
    plot_st <- st_list_NEA_cl[[y]][[m]]
    
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
  
  #crop to extent of BPNS
  pr_geo_y <- crop(pr_geo_y, bpns_shp)
  pr_geo_y <- mask(pr_geo_y, bpns_shp)
  
  df <- values(pr_geo_y) %>%
    as.data.frame() %>%
    na.omit() %>%
    gather(key = month, value = "value") %>%
    mutate(month = str_remove(month,"X"))
  
  df_out <- cbind(df_out, df[,2])
  names(df_out)[y+1] <- as.character(c(2000:2020)[y])
  
  print(c(2000:2020)[y])
}

df_evol <- df_out %>%
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

ggsave(paste0("3.MODEL_OUTPUT/ADULTS/adults_average_monthly_variation_HSI.png"), width = 10, height = 4)
