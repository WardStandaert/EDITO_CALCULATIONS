#install and load packages
pckgs <- c("raster","sp","proj4","ncdf4","car","mgcv","dismo","rJava","ENMeval",
           "boot","gstat","mgcv", "ggplot2", "tidyr", "dynamicSDM", "stringr", 
           "mapdata", "base","tidync", "sf", "mapview", "tictoc", "ape", "spdep", 
           "spThin", "StatMatch", "CoordinateCleaner", "maxnet", "rasterVis",
           "tibble", "purrr", "tidyverse", "arrow")

installed_packages <- pckgs %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  install.packages(pckgs[!installed_packages])
}

invisible(lapply(pckgs, library, character.only = TRUE))
rm(pckgs)

#Rarr package needs separate install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rarr")
library("Rarr")

##########  PART 1 - train model  ########## 

setwd("/home/onyxia/work/EDITO_CALCULATIONS/")
#1. Collect & organize occurrence data -----

# the file to process
acf <- S3FileSystem$create(
  anonymous = T,
  scheme = "https",
  endpoint_override = "s3.waw3-1.cloudferro.com"
)

eurobis <- arrow::open_dataset(acf$path("emodnet/biology/eurobis_occurence_data/eurobisgeoparquet/eurobis_no_partition_sorted.parquet" ))
df_occs <- eurobis |> 
  filter(aphiaidaccepted==126417, datasetid==4423,
         longitude > -12, longitude < 10,
         latitude > 48, latitude < 62,
         observationdate >= as.POSIXct("2000-01-01"),
         observationdate <= as.POSIXct("2020-12-31")) |> 
  collect() 

glimpse(df_occs)

df_occs <- df_occs %>% 
  select(Latitude=latitude,
         Longitude=longitude, 
         Time=observationdate) %>%
  mutate(year = year(Time),
         month = month(Time))

mapview(df_occs %>% dplyr::select(Longitude) %>% pull,
        df_occs %>% dplyr::select(Latitude) %>% pull,
        crs = "epsg:4326")

table(df_occs$month)

# remove outliers
# df_occs2 <- CoordinateCleaner::cc_outl(x = as.data.frame(df_occs), lon = "lon", lat = "lat", 
#                                        method = "quantile", mltpl = 1.5, verbose = TRUE)

#Removed 0 records.


#2. Remove duplicates ----
nrow(df_occs)
# 7902
df_occs <- df_occs %>%
  distinct(year, month, Longitude, Latitude, .keep_all = TRUE)
nrow(df_occs)
# 5347
# so we removed 2555 duplicated records


#3. Reduce sampling bias ----
# Thin towards distance of 10 NM or 18.52 km
# This is the distance between valid hauls in ICES trawl surveys

df_occs_thinned <- df_occs[0,] %>% select(-Time)
for (y in 1:length(2000:2020)) {
  for (m in 1:length(1:12)) {
    tmp_df <- df_occs %>% filter(year == c(2000:2020)[y],
                                 month == m)
    tmp_df_thinned <- spThin::thin(tmp_df %>% mutate(species = "Atlantic herring"),
                                   lat.col = "Latitude", long.col = "Longitude",
                                   spec.col = "species", thin.par = 18.52, reps = 1,
                                   write.files = FALSE, locs.thinned.list.return = TRUE,
                                   verbose = FALSE)[[1]]
    
    tmp_df_thinned <- tmp_df_thinned %>%
      mutate(year = c(2000:2020)[y],
             month = m)
    
    df_occs_thinned <- rbind(df_occs_thinned, tmp_df_thinned)
  }
  print(paste0(c(2000:2020)[y], " done"))
}


nrow(df_occs_thinned)
# 3995

mapview(df_occs_thinned %>% dplyr::select(Longitude) %>% pull,
        df_occs_thinned %>% dplyr::select(Latitude) %>% pull,
        crs = "epsg:4326")

# save(df_occs_thinned, file = "data/df_thinned.Rdata")
load("data/df_thinned.Rdata")

#4. Create background points ----
# lets take a spatial buffer of 100 km around the occurrence points for this example 
# 
# abs <- spatiotemp_pseudoabs(spatial.method = "buffer", temporal.method = "random",
#                             occ.data = df_occs_thinned %>% mutate(x = Longitude, y = Latitude),
#                             temporal.ext = c("2000-01-01", "2020-12-31"), spatial.buffer = 100000,
#                             n.pseudoabs = 10000)
# 
# glimpse(abs)
# 
# #limit temporal values of background points to months where occurrences of larvae are present
# set.seed(123)
# abs <- abs %>%
#   mutate(month = sample(unique(df_occs_thinned$month), nrow(abs), replace = TRUE)) %>%
#   mutate(Longitude = x, Latitude = y) %>%
#   select(-day, -x, -y)
# glimpse(abs)
# 
# save(abs, file = "zarr_extraction/absences_save.Rdata")
load("zarr_extraction/absences_save.Rdata")

#combine occurrences and background points
df_occ_bg <- rbind(df_occs_thinned %>% mutate(presence = 1), 
                   abs %>% mutate(presence = 0))
glimpse(df_occ_bg)
save(df_occ_bg, file = "zarr_extraction/pres_abs_save.Rdata")

#5. Sample environmental variables at occurrences and background points ----
source("zarr_extraction/editoTools.R")
options("outputdebug"=c('L','M'))
load(file = "zarr_extraction/editostacv2.par")

#the requested timestep resolution of the dataset in milliseconds
#in this case we work with monthly data (1 month = 30.436875*24*3600*1000 = 2629746000 milliseconds)
timeSteps=c(2629746000)

##TODO: ADD ELEVATION, WINDFARMS AND SEABED SUBSTRATE ENERGY ------
parameters = list("thetao" = c("par" = "thetao", "fun" = "mean", "buffer" = "18520"),
                  "so" = c("par" = "so", "fun" = "mean", "buffer" = "18520"),
                  "zooc" = c("par" = "zooc", "fun" = "mean", "buffer" = "18520"),
                  "phyc" = c("par" = "phyc", "fun" = "mean", "buffer" = "18520"))

#check if they are all available in the data lake
for (parameter in parameters) {
  param = ifelse(length(parameter) == 1, parameter, parameter["par"])
  if(! param %in% unique(EDITOSTAC$par)) dbl("Unknown parameter ", param)
}

#extract function (requires POSIXct Time column)
df_occ_bg_env = enhanceDF(inputPoints = df_occ_bg %>% 
                            mutate(Time = as.POSIXct(paste(year,month,1,sep = "-"))),
                          requestedParameters = parameters,
                          requestedTimeSteps = timeSteps,
                          stacCatalogue = EDITOSTAC,
                          verbose="on")

glimpse(df_occ_bg_env)

df_occ_bg_env <- df_occ_bg_env %>% select(Longitude, Latitude, year, month,
                                          thetao, so, zooc, phyc)

# save(df_occ_bg_env, file = "zarr_extraction/pres_abs_env_save.Rdata")
load(file = "zarr_extraction/pres_abs_env_save.Rdata")


# remove observations where no environmental values were present
df_occ_bg_env <- drop_na(df_occ_bg_env)
table(df_occ_bg_env$presence)
#    0    1 
# 7617 3910 


#6 ENMeval model creation  ----

glimpse(df_occ_bg_env)

## TODO remove this part ----
substr_lvl <- tibble(sub_char = c("Fine mud", "Sand", "Muddy sand", "Mixed sediment",
                                  "Coarse substrate","Sandy mud or Muddy sand", "Seabed",
                                  "Rock or other hard substrata","Sandy mud", "Sandy mud or Muddy sand ",
                                  "Sediment","Fine mud or Sandy mud or Muddy sand"),
                     seabed_substrate = c(1:12))
energy_lvl <- tibble(ene_char = c("High energy", "Moderate energy", "Low energy", "No energy information"),
                     seabed_energy = c(1:4))

df_occ_bg_env <- df_occ_bg_env %>%
  left_join(energy_lvl, by = "seabed_energy") %>%
  left_join(substr_lvl, by = "seabed_substrate") %>%
  dplyr::select(-seabed_energy, -seabed_substrate)


#if any column exists with only one unique value --> delete it
if(any(apply(rbind(tmp_pr_predv, tmp_bg_predv), 2, function(x) length(unique(x)) == 1))) {
  tmp_ind <- which(apply(rbind(tmp_pr_predv, tmp_bg_predv), 2, function(x) length(unique(x)) == 1))
  tmp_ind_v <- names(tmp_ind)
  tmp_pr_predv <- tmp_pr_predv[-tmp_ind]
  tmp_bg_predv <- tmp_bg_predv[-tmp_ind]
  v_list <- v_list[-which(v_list == tmp_ind_v)]
  print(paste0("removed ", names(tmp_ind), " as a variable for ", species_name$simple))
}

# input for ENMevaluate requires lon & lat as first covariates

#test different model settings using ENMeval
model_fit <- ENMeval::ENMevaluate(occs = df_occ_bg_env %>% 
                                    filter(presence == 1) %>% 
                                    select(Longitude, Latitude, depth, Phyto, SSS, SST, windfarms, ZooPl, ene_char, sub_char) %>%
                                    mutate(windfarms = as.factor(windfarms)),
                                  bg = df_occ_bg_env %>% 
                                    filter(presence == 0) %>% 
                                    select(Longitude, Latitude, depth, Phyto, SSS, SST, windfarms, ZooPl, ene_char, sub_char) %>%
                                    mutate(windfarms = as.factor(windfarms)),
                                  tune.args = list(fc = c("L","LQ","LQH"),
                                                   rm = c(1,2,4,8,32)),
                                  algorithm = "maxnet",
                                  partitions = "randomkfold",
                                  categoricals = c("sub_char", "ene_char", "windfarms"),
                                  doClamp = TRUE,
                                  parallel = TRUE)




print(model_fit %>% eval.results() %>% filter(delta.AICc == 0))

save(model_fit, file = "model_fit_save.Rdata")
load("model_fit_save.Rdata")


#7 model evaluation AUC & TSS ----
AUC_maxent <- list()
TSS_maxent <- list()

Presences <- pr_predv
Background <- bg_predv

if(all(c("seabed_energy","seabed_substrate") %in% names(pr_predv))) {
  Presences <- pr_predv %>%
    left_join(energy_lvl, by = "seabed_energy") %>%
    left_join(substr_lvl, by = "seabed_substrate") %>%
    dplyr::select(-seabed_energy, -seabed_substrate)
  
  Background <- bg_predv %>%
    left_join(energy_lvl, by = "seabed_energy") %>%
    left_join(substr_lvl, by = "seabed_substrate") %>%
    dplyr::select(-seabed_energy, -seabed_substrate)
}
else {
  Presences <- pr_predv
  Background <- bg_predv
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
  eval_res <- eval_res_list_sp
  opt.aicc <- eval.results(eval_res_list_sp) %>% filter(delta.AICc == 0)
  mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]
  
  EvalBgRes  <- predict(mod, EvalBg, type = "cloglog")
  EvalPresRes  <- predict(mod, EvalPres, type = "cloglog")
  
  tmp_AUC[[i]]<-evaluate(c(EvalPresRes), c(EvalBgRes))
  umbral_maxent<-threshold( tmp_AUC[[i]])
  tmp_TSS[[i]]<-evaluate(p=c(EvalPresRes), a=c(EvalBgRes), tr=umbral_maxent$spec_sens)
  
  toc()
}

AUC_maxent <- tmp_AUC
TSS_maxent <- tmp_TSS
names(AUC_maxent) <- names(TSS_maxent) <- species_name$simple


AUC_maxent_vect   <- data.frame(matrix(nrow = 10, ncol = 2))
TPR_maxent_vect   <- data.frame(matrix(nrow = 10, ncol = 2))
TNR_maxent_vect   <- data.frame(matrix(nrow = 10, ncol = 2))
TSS_maxent_vect   <- data.frame(matrix(nrow = 10, ncol = 2))

for (s in 1:1) {
  AUC_maxent_vect[,s]  <- sapply(AUC_maxent,function(x){slot(x,'auc')})
  TPR_maxent_vect[,s]  <- sapply(TSS_maxent,function(x){slot(x,'TPR')})
  TNR_maxent_vect[,s]  <- sapply(TSS_maxent,function(x){slot(x,'TNR')})
  TSS_maxent_vect[,s]  <- TPR_maxent_vect[,s] + TNR_maxent_vect[,s] - 1
  
  colnames(AUC_maxent_vect) <- colnames(TPR_maxent_vect) <-
    colnames(TNR_maxent_vect) <-  colnames(TSS_maxent_vect) <-
    species_name$simple
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
sum_all[match(sum_all$species, species_name$simple),]

#save results
write.csv(sum_all, "3.MODEL_OUTPUT/LARVAE/validation_metrics.csv")


##########  PART 2 - project model  ########## 

#9. Spatial predictions ----
for (m in 1:12) {
  pr_geo_y <- stack()
  for (y in 1:length(2000:2020)) {
    #get model
    eval_res <- eval_res_list_sp[[1]]
    opt.aicc <- eval.results(eval_res_list_sp[[1]]) %>% filter(delta.AICc == 0)
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
  
  writeRaster(av_NEA, paste0("3.MODEL_OUTPUT/LARVAE/RASTERS/larval_herring_average_HSI_", m, ".tif"), overwrite = T)
  writeRaster(sd_NEA, paste0("3.MODEL_OUTPUT/LARVAE/RASTERS/larval_herring_sd_HSI_", m, ".tif"), overwrite = T)
  
  
  if (m == 10) {
    # variability BPNS
    bpns_shp <- st_read("1.DOWNLOAD/shapefiles/BPNS/eez.shp", quiet = TRUE)
    
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
    
    ggsave(paste0("3.MODEL_OUTPUT/LARVAE/larvae_yearly_variation_", m,".png"), width = 10, height = 4)
  }
  print(m)
}

