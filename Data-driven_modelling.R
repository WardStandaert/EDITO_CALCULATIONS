pckgs <- c("raster","sp","proj4","ncdf4","car","mgcv","dismo","rJava","ENMeval",
           "boot","gstat","mgcv", "ggplot2", "tidyr", "dynamicSDM", "stringr", 
           "mapdata", "base","tidync", "sf", "mapview", "tictoc", "ape", "spdep", 
           "spThin", "StatMatch", "CoordinateCleaner", "maxnet", "rasterVis",
           "tibble", "purrr", "tidyverse")
pckgs[which(lapply(pckgs, require, character.only = TRUE) == FALSE)]
rm(pckgs)

setwd("/home/onyxia/work/marbcubes/")
#1. Read & organize data -----
df_fish <- read.csv("data/raw_data/bcube_fcomm.csv")
df_fish <- df_fish %>%
  dplyr::mutate(Longitude = as.numeric(map_chr(str_split(cellCode, "_"),2)),
                Latitude = as.numeric(map_chr(str_split(cellCode, "_"),1)),
                year = as.numeric(map_chr(str_split(yearMonth, "-"),1)),
                month = as.numeric(map_chr(str_split(yearMonth, "-"),2)),
                Time = as.POSIXct(paste(1, month, year, sep = "-"), format = "%d-%m-%Y")) %>%
  dplyr::filter(year >= 2010 & year <= 2020)

glimpse(df_fish)

mapview(df_fish %>% dplyr::select(Longitude) %>% pull,
        df_fish %>% dplyr::select(Latitude) %>% pull,
        crs = "epsg:4326")

#filter out the spatial extent of the fish' native area
df_fish <- df_fish %>%
  filter(Longitude >= 30,
         Longitude <= 160,
         Latitude >= -35,
         Latitude <= 32) %>%
  mutate(species = "cornetfish")

mapview(df_fish %>% dplyr::select(Longitude) %>% pull,
        df_fish %>% dplyr::select(Latitude) %>% pull,
        crs = "epsg:4326")

table(df_fish$month)

#2. Remove duplicates ----
df_fish <- df_fish %>%
  distinct(year, month, Longitude, Latitude, .keep_all = TRUE)


#3. create pseudo-absences ----
abs <- spatiotemp_pseudoabs(spatial.method = "buffer", temporal.method = "random",
                            occ.data = fish_df2 %>% mutate(x = Longitude, y = Latitude),
                            temporal.ext = c("2010-01-01", "2020-12-31"),
                            spatial.buffer = 10000, n.pseudoabs = 1000)

df_p <- df_fish %>% 
  select(Longitude, Latitude, year, month, Time) %>%
  mutate(pa = 1)
df_a <- abs %>% rename(Longitude = x, Latitude = y) %>%
  mutate(Time = as.POSIXct(paste(year,month,day, sep = "-")),
         pa = 0) %>%
  select(-day)


#4. Sample environmental variables with biological data ----
source("src/data-driven_model_fcomm/editoTools.R")
options("outputdebug"=c('L','M'))

#the cached stacCatalog is called 'EDITOSTAC'
load(file = "src/data-driven_model_fcomm/editostacv2.par")

#the requested timestep resolution of the dataset in milliseconds
#in this case we work with monthly data (1 month = 30.436875*24*3600*1000 = 2629746000 milliseconds)
timeSteps=c(2629746000)

parameters_pres = list("elevation" = c("par" = "elevation",
                                  "fun" = "mean",
                                  "buffer" = "10000"),
                  "thetao"= c("par" = "thetao",
                              "fun" = "mean",
                              "buffer" = "10000"),
                  "so"= c("par" = "so",
                          "fun" = "mean",
                          "buffer" = "10000"),
                  "chl"= c("par" = "chl",
                              "fun" = "mean",
                              "buffer" = "10000"))

parameters_abs = list("elevation",
                       "thetao",
                       "so",
                       "chl")

#check if they all exist
for ( parameter in parameters_pres) {
  param = ifelse(length(parameter) == 1, parameter, parameter["par"])
  if(! param %in% unique(EDITOSTAC$par) )
  { dbl("Unknown parameter ", param)
  }
}

# Extract data ----
#add verbose= anything to get additional info on the positions ( par_x, par_y, par_z ) and time (par_t) found in the zarr files

enhanced_DF_pres = enhanceDF(inputPoints = df_p,
                             requestedParameters = parameters_pres,
                             requestedTimeSteps = timeSteps,
                             stacCatalogue = EDITOSTAC,
                             verbose="on")

glimpse(enhanced_DF_pres)

enhanced_DF_abs = enhanceDF(inputPoints = df_a,
                             requestedParameters = parameters_abs,
                             requestedTimeSteps = timeSteps,
                             stacCatalogue = EDITOSTAC,
                             verbose="on")

glimpse(enhanced_DF_abs)

fish_df2 <- rbind(enhanced_DF_pres %>% 
                    select(Longitude, Latitude, year, month,
                           Time, SST = thetao, SSS = so, chl) %>%
                    mutate(species = "fish", pa = 1),
                  enhanced_DF_abs %>% 
                    select(Longitude, Latitude, year, month,
                           Time, SST = thetao, SSS = so, chl)%>%
                    mutate(species = "fish", pa = 0))

nrow(fish_df2)
fish_df2 <- na.omit(fish_df2)
nrow(fish_df2)

#3. Sampling bias - filtering -----

### sampling bias: spatial filtering ----
#already thinned to minimum distance of 10 km


# Variables selected through literature review
v_list <- c("SST", "SSS", "chl")

#5. extract environmental values for presence absences ----




#6 ENMeval model creation  ----
#convert month to character bc algorithm cant handle factors that look like numeric

eval_res_list_sp <- list()
tmp_pr_predv <- list()
tmp_bg_predv <- list()

if(all(c("seabed_energy","seabed_substrate") %in% names(pr_predv))) {
  tmp_pr_predv <- pr_predv %>%
    left_join(energy_lvl, by = "seabed_energy") %>%
    left_join(substr_lvl, by = "seabed_substrate") %>%
    dplyr::select(-seabed_energy, -seabed_substrate)

  tmp_bg_predv <- bg_predv %>%
    left_join(energy_lvl, by = "seabed_energy") %>%
    left_join(substr_lvl, by = "seabed_substrate") %>%
    dplyr::select(-seabed_energy, -seabed_substrate)
  cat <- c("ene_char","sub_char")
} else {
  tmp_pr_predv <- pr_predv
  tmp_bg_predv <- bg_predv
  cat <- NULL
}

#if any column exists with only one unique value --> delete it
if(any(apply(rbind(tmp_pr_predv, tmp_bg_predv), 2, function(x) length(unique(x)) == 1))) {
  tmp_ind <- which(apply(rbind(tmp_pr_predv, tmp_bg_predv), 2, function(x) length(unique(x)) == 1))
  tmp_ind_v <- names(tmp_ind)
  tmp_pr_predv <- tmp_pr_predv[-tmp_ind]
  tmp_bg_predv <- tmp_bg_predv[-tmp_ind]
  v_list <- v_list[-which(v_list == tmp_ind_v)]
  print(paste0("removed ", names(tmp_ind), " as a variable for ", species_name$simple))
}

tmp_pr_predv <- data.frame(pr_coord %>% dplyr::select(lon, lat), tmp_pr_predv) #input for ENMevaluate requires lon & lat as first covariates
tmp_bg_predv <- data.frame(bg_coord %>% dplyr::select(lon, lat), tmp_bg_predv) #input for ENMevaluate requires lon & lat as first covariates

#test different model settings using ENMeval
eval_res_list_sp <- ENMeval::ENMevaluate(occs = tmp_pr_predv,
                                         bg = tmp_bg_predv,
                                         tune.args = list(fc = c("L","LQ","LQH"),
                                                          rm = c(1,2,4,8,32)),
                                         algorithm = "maxnet",
                                         partitions = "randomkfold",
                                         categoricals = cat,
                                         doClamp = TRUE,
                                         parallel = TRUE)



print(eval_res_list_sp %>% eval.results() %>% filter(delta.AICc == 0))


#7 model evaluation AUC & TSS ----
set.seed(123)

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


#8. variable importance ----

# make big tibble with all possible values per variable
# only run once for either adults or larvae

# tib_all_v <- tibble(a = rep(1,3846528)) %>% select()
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
# tib_all_v$Year <- rep(2000:2020, times = 1, each = ncell(st_list_NEA2[[1]][[1]])*12)
# colnames(tib_all_v)[which(colnames(tib_all_v) %in% c("ene_char", "sub_char"))] <- c("seabed_energy", "seabed_substrate")
#
# tib_all_v <- tib_all_v %>%
#   left_join(substr_lvl %>% mutate(seabed_substrate = as.numeric(seabed_substrate)), by = "seabed_substrate") %>%
#   left_join(energy_lvl %>% mutate(seabed_energy = as.numeric(seabed_energy)), by = "seabed_energy") %>%
#   dplyr::select(-seabed_substrate, -seabed_energy)
#
#
# head(tib_all_v)
#
# save(tib_all_v, file = "SAVE/tib_all_v.Rdata")
load("SAVE/tib_all_v.Rdata")

cor_maxnet <- list()
#1. get best model (based on AICc criterion)
eval_res <- eval_res_list_sp
res <- eval.results(eval_res)
opt.aicc <- res %>% filter(delta.AICc == 0)
mod.seq <- eval.models(eval_res)[[opt.aicc$tune.args]]

#2. make original prediction
df_pred <- rbind(pr_predv, bg_predv)
if(all(c("seabed_energy", "seabed_substrate") %in% names(df_pred))) {
  df_pred <- df_pred %>%
    left_join(substr_lvl, by = "seabed_substrate") %>%
    left_join(energy_lvl, by = "seabed_energy") %>%
    select(-seabed_substrate, -seabed_energy)
}

prediction1 <- predict(mod.seq, df_pred,
                       se.fit=TRUE, type = "cloglog",
                       factors = list(seabed_energy = factor(energy_lvl$seabed_energy,
                                                             labels = energy_lvl$seabed_energy,
                                                             levels = energy_lvl$seabed_energy),
                                      seabed_substrate = factor(substr_lvl$seabed_substrate,
                                                                labels = substr_lvl$seabed_substrate,
                                                                levels = substr_lvl$seabed_substrate)))

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
    prediction2 <- predict(mod.seq, df_pred_res, se.fit=TRUE, type = "cloglog")
    cor_new_max <- cor(na.omit(prediction1), na.omit(prediction2))

    #calculate correlation between the two predictions
    cor_df_max[i, v] <- cor_new_max
  }
  toc()
}
#store in list
cor_maxnet <- sapply(cor_df_max, FUN = function(x) 1-mean(x, na.rm=T))


cor_maxnet_sp <- cor_maxnet

#convert to percentages
var_imp_lv <- cor_maxnet_sp[[1]]
var_imp_lv <- round(var_imp_lv/sum(var_imp_lv),2)
var_imp_lv
sum(var_imp_lv)

write.csv(var_imp_lv, "3.MODEL_OUTPUT/LARVAE/variable_importance.csv")

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

