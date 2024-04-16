#install and load packages
pckgs <- c("raster","sp","proj4","ncdf4","car","mgcv","dismo","rJava","ENMeval",
           "boot","gstat","mgcv", "ggplot2", "tidyr", "dynamicSDM", "stringr", 
           "mapdata", "base","tidync", "sf", "mapview", "tictoc", "ape", "spdep", 
           "spThin", "StatMatch", "CoordinateCleaner", "maxnet", "rasterVis",
           "tibble", "purrr", "tidyverse", "arrow", "tidyverse")

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
load(file = "pres_abs_env_save.Rdata")


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
  dplyr::select(-seabed_energy, -seabed_substrate) %>%
  mutate(windfarms = as.factor(windfarms))


# input for ENMevaluate requires lon & lat as first covariates

#test different model settings using ENMeval
model_fit <- ENMeval::ENMevaluate(occs = df_occ_bg_env %>%                                    
                                    filter(presence == 1) %>% 
                                    dplyr::select(Longitude, Latitude, depth, Phyto, SSS, SST, windfarms, ZooPl, ene_char, sub_char),
                                  bg = df_occ_bg_env %>% 
                                    filter(presence == 0) %>% 
                                    dplyr::select(Longitude, Latitude, depth, Phyto, SSS, SST, windfarms, ZooPl, ene_char, sub_char),
                                  tune.args = list(fc = c("L","LQ","LQH"),
                                                   rm = c(1,2,4,8,32)),
                                  algorithm = "maxnet",
                                  partitions = "randomkfold",
                                  categoricals = c("sub_char", "ene_char", "windfarms"),
                                  doClamp = TRUE,
                                  parallel = TRUE)




print(model_fit %>% eval.results() %>% filter(delta.AICc == 0))

# save(model_fit, file = "model_fit_save.Rdata")
load("model_fit_save.Rdata")


#7 model evaluation AUC & TSS ----
AUC_maxent <- list()
TSS_maxent <- list()

Presences <- df_occ_bg_env %>%                                    
  filter(presence == 1) %>% 
  dplyr::select(Longitude, Latitude, depth, Phyto, SSS, SST, windfarms, ZooPl, ene_char, sub_char)
Background <- df_occ_bg_env %>%                                    
  filter(presence == 0) %>% 
  dplyr::select(Longitude, Latitude, depth, Phyto, SSS, SST, windfarms, ZooPl, ene_char, sub_char)

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
  eval_res <- model_fit
  opt.aicc <- eval.results(model_fit) %>% filter(delta.AICc == 0)
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

AUC_maxent_vect  <- sapply(AUC_maxent,function(x){slot(x,'auc')})
TPR_maxent_vect  <- sapply(TSS_maxent,function(x){slot(x,'TPR')})
TNR_maxent_vect  <- sapply(TSS_maxent,function(x){slot(x,'TNR')})
TSS_maxent_vect  <- TPR_maxent_vect + TNR_maxent_vect - 1

stat_df <- data.frame(AUC = c(mean(AUC_maxent_vect), sd(AUC_maxent_vect)),
                      TSS = c(mean(TSS_maxent_vect), sd(TSS_maxent_vect)),
                      TPR = c(mean(TPR_maxent_vect), sd(TPR_maxent_vect)),
                      TNR = c(mean(TNR_maxent_vect), sd(TNR_maxent_vect)))

stat_df <- data.frame(statistic = c("AUC", "TSS", "TPR", "TNR"),
                      mean = c(mean(AUC_maxent_vect), mean(TSS_maxent_vect), mean(TPR_maxent_vect), mean(TNR_maxent_vect)),
                      sd = c(sd(TSS_maxent_vect), sd(TSS_maxent_vect), sd(TPR_maxent_vect), sd(TNR_maxent_vect)))
  
stat_df

##########  PART 2 - project model  ########## 

#9. get raster slices to project on ----
variables <- c("depth", "SST", "SSS", "Phyto",
               "ZooPl", "windfarms", "seabed_substrate", 
               "seabed_energy")

files <- tibble(name = list.files("data-raw", pattern = ".tif|.grd", full.names = TRUE),
                parameter = str_extract(name, pattern = paste0(variables, collapse = "|")),
                month = str_extract(name, pattern = "\\d{2}|\\d"))

#get model
eval_res <- model_fit
opt.aicc <- eval.results(model_fit) %>% filter(delta.AICc == 0)
mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args)]]

#10. Spatial predictions ----

predictions <- stack()
for (m in unique(df_occ_bg_env$month)) {
  
  plot_st <- stack(files %>% 
                     filter(month == m | is.na(month)) %>% 
                     dplyr::select(name) %>% 
                     pull())
  #get model
  
  names(plot_st) <- str_extract(names(plot_st), pattern = paste0(variables, collapse = "|"))
  names(plot_st)[which(names(plot_st) == "seabed_energy")] <- "ene_char"
  names(plot_st)[which(names(plot_st) == "seabed_substrate")] <- "sub_char"
  
  pred_m <- predict(plot_st, mod, clamp=T, type="cloglog",
                         factors = list(ene_char = factor(energy_lvl$seabed_energy,
                                                          labels = energy_lvl$ene_char,
                                                          levels = energy_lvl$seabed_energy),
                                        sub_char = factor(substr_lvl$seabed_substrate,
                                                          labels = substr_lvl$sub_char,
                                                          levels = substr_lvl$seabed_substrate)))
  # crop to extent of study area
  pred_m <- crop(pred_m, extent(-12, 10, 48, 62))
  
  names(pred_m) <- paste0("prediction_", month.abb[m])

  predictions <- stack(predictions, pred_m)
  
  print(m)
}

plot(predictions)


#11. Response curves ----
length(variables)
vs <- c("depth", "SST", "SSS", "Phyto", "ZooPl", "windfarms", "sub_char", "ene_char")

par(mfrow = c(3,3))
for (v in vs) {
  response.plot(mod, v, type = "cloglog", 
                ylab = "Probability of occurrence",
                plot = T)
}


#12. Add variable importance? ----
