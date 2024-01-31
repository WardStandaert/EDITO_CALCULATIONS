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

#1.general setup ----
setwd("DATA")
load("SAVE/final_st_list.Rdata")
load("SAVE/bio_dat_larvae.Rdata")
load("SAVE/presences_after_spatial_filtering.Rdata")
rm(presences_list)

all_vars <- str_remove_all(names(st_list_NEA_cl[[1]][[1]]) , "_\\d{4}_\\d{2}|_\\d{4}_\\d")

species_sp <- tibble(scientific = c("Clupea harengus", "Scomber scombrus"),
                  simple = c("herring", "mackerel"))


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


#2. specific setup ----
runX <- 1
txt <- "number of pres lots of variations"
varying <- "npres"
n_pres_pts <- rep(seq(from = 50, to = 100, by = 50),1)  #numeric vector of length 4

v_list <- list()
v_list$herring <- c("SST", "SSS", "windfarms", "Phyto", "ZooPl", "seabed_energy", "seabed_substrate", "depth")
v_list$mackerel <- c("SST", "SSS", "windfarms", "max_SSV", "ZooPl", "depth")

# v_list2 <- list()
# v_list2$herring <- c("Year","Month","SST", "SSS", "windfarms", "Phyto", "ZooPl", "seabed_energy", "seabed_substrate", "depth")
# v_list2$mackerel <- c("Year","Month","SST", "SSS", "windfarms", "max_SSV", "EuphD", "ZooPl", "depth")
# v_list_specific <- list(v_list, v_list2)

fig_out_name <- paste0("_response_curve_", varying, "_")

dir.create(paste0("SAVE/test_runs_larvae/run_", runX), showWarnings = TRUE)
write.table(txt, paste0("SAVE/test_runs_larvae/run_", runX, "/config.txt"))

df_out <- data.frame(species = rep(species_sp$simple, length(n_pres_pts)),
                     n_pres_pts = rep(n_pres_pts, each = nrow(species_sp)))
# df_out <- data.frame(species = rep(species$simple, length(n_pres_pts)),
                     # n_pres_pts = rep(n_pres_pts, each = nrow(species)),
                     # vars = rep(c("vlist1","vlist2"), each = nrow(species)))

df_out$SACp_pval <- df_out$SACp <- df_out$AUC <- df_out$rm <- df_out$fc <- df_out$n_bg_pts <- NA

for (n in 1:length(n_pres_pts)){
  # v_list <- v_list_specific[[n]]
  tmp_n_pres <- n_pres_pts[n]
  
  #2.1 environmental sampling ----
  red_dat <- list()
  for (s in 1:nrow(species_sp)) {
    name <- species_sp$simple[s]

    tmp_full_dat <- left_join(presences_list_sp[[which(names(presences_list_sp) == species_sp$simple[s])]] %>%
                                mutate(lon = Longitude, lat = Latitude), full_dat_sp, by = c("lon","lat")) %>%
      filter(scientificname == species_sp$scientific[s])

    df_pr <- tmp_full_dat[0,]
    df_rem <- tmp_full_dat

    #remove variables that have only one unique value
    tmp_df <- tmp_full_dat %>%
      select(any_of(all_vars)) %>%
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
    while (nrow(df_pr) < tmp_n_pres & nrow(df_rem) > 0) {
      mal_dis <- mahalanobis.dist(df_pr[,12:16],df_rem[,12:16])
      mal_dis_tot <- colSums(mal_dis)
      id <- as.vector(which(mal_dis_tot == max(mal_dis_tot), useNames = FALSE))
      if (length(id) > 1) {id <- sample(id,1)}
      df_pr <- rbind(df_pr,df_rem[id,])
      df_rem <- df_rem[-id,]
    }
    red_dat[[s]] <- df_pr
  }
  reduced_dat <- bind_rows(red_dat)
  
  #2.1 without environmental sampling ----
  # red_dat <- list()
  # for (s in 1:nrow(species_sp)) {
  #   name <- species_sp$simple[s]
  #   
  #   tmp_full_dat <- left_join(presences_list_sp[[which(names(presences_list_sp) == species_sp$simple[s])]] %>%
  #                               mutate(lon = Longitude, lat = Latitude), full_dat_sp, by = c("lon","lat")) %>%
  #     filter(scientificname == species_sp$scientific[s])
  #   
  #   df_pr <- tmp_full_dat
  #   red_dat[[s]] <- df_pr
  # }
  # reduced_dat <- bind_rows(red_dat)
  
  print("1.filtering done")
  
  #2.2 presence & background points ----
  pr_pa <- list()
  pr_coord <- list()
  pr_predv <- list()
  
  for (s in 1:nrow(species_sp)) {
    tmp_pres  <- reduced_dat %>%
      filter(scientificname == species_sp$scientific[s])
    
    pr_pa[[s]] <- rep(1, nrow(tmp_pres))
    pr_coord[[s]] <- tmp_pres %>% select(lon, lat, Year, Month)
    pr_predv[[s]] <- tmp_pres %>% select(all_of(v_list[[which(names(v_list) == species_sp$simple[s])]]))
    
    names(pr_pa)[s] <- names(pr_coord)[s] <- names(pr_predv)[s] <- species_sp$simple[s]
  }
  
  #background points
  
  bg_pa <- list()
  bg_coord <- list()
  bg_predv <- list()
  
  for (s in 1:nrow(species_sp)) {
    #create background points within ices regions that species_sp occurs in
    tmp_dat <- reduced_dat_sp %>%
      filter(scientificname == species_sp$scientific[s])
    
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
                                       data.frame(tmp_coords_dat$lon, tmp_coords_dat$lat), method = "simple")
        colnames(tmp_extract) <- colnames(tmp_extract) %>%
          str_remove_all("_\\d{4}_\\d{2}|_\\d{4}_\\d")
        tmp_df <- data.frame(tmp_coords_dat, tmp_extract)
        tmp_bckgr <- rbind(tmp_bckgr, tmp_df)
      }
    }
    
    #outputs
    bg_pa[[s]] <- rep(0, nrow(tmp_bckgr))
    
    bg_coord[[s]] <- tmp_bckgr %>% select(lon, lat, Year, Month)
    
    bg_predv[[s]] <- tmp_bckgr %>% 
      select(all_of(v_list[[which(names(v_list) == species_sp$simple[s])]])) #variable selection
    
    names(bg_pa)[s] <- names(bg_coord)[s] <- names(bg_predv)[s] <- species_sp$simple[s]
  }
  
  #remove all tmp files
  rm(list = ls()[which(str_detect(ls(), pattern = "^tmp_"))])
  
  
  print("2.presence & background points created & sampled")
  
  #2.3. ENMevaluate model creation ----
  eval_res_list <- list()
  tmp_pr_predv <- list()
  tmp_bg_predv <- list()
  for (s in 1:nrow(species_sp)) {
    if(all(c("seabed_energy","seabed_substrate","Month") %in% names(pr_predv[[s]]))) {
      tmp_pr_predv[[s]] <- pr_predv[[s]] %>%
        left_join(month_lvl, by = "Month") %>%
        left_join(energy_lvl, by = "seabed_energy") %>%
        left_join(substr_lvl, by = "seabed_substrate") %>%
        select(-Month, -seabed_energy, -seabed_substrate)
      
      tmp_bg_predv[[s]] <- bg_predv[[s]] %>%
        left_join(month_lvl, by = "Month") %>%
        left_join(energy_lvl, by = "seabed_energy") %>%
        left_join(substr_lvl, by = "seabed_substrate") %>%
        select(-Month, -seabed_energy, -seabed_substrate)
      cat <- c("mon_char", "ene_char","sub_char")
    }
    else if(all(c("seabed_energy","seabed_substrate") %in% names(pr_predv[[s]]))) {
      tmp_pr_predv[[s]] <- pr_predv[[s]] %>%
        left_join(energy_lvl, by = "seabed_energy") %>%
        left_join(substr_lvl, by = "seabed_substrate") %>%
        select(-seabed_energy, -seabed_substrate)
      
      tmp_bg_predv[[s]] <- bg_predv[[s]] %>%
        left_join(energy_lvl, by = "seabed_energy") %>%
        left_join(substr_lvl, by = "seabed_substrate") %>%
        select(-seabed_energy, -seabed_substrate)
      cat <- c("ene_char","sub_char")
    }
    else if("Month" %in% names(pr_predv[[s]])) {
      tmp_pr_predv[[s]] <- pr_predv[[s]] %>%
        left_join(month_lvl, by = "Month") %>%
        select(-Month)
      
      tmp_bg_predv[[s]] <- bg_predv[[s]] %>%
        left_join(month_lvl, by = "Month") %>%
        select(-Month)
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
      print(paste0("removed ", names(tmp_ind), " as a variable for ", species_sp$simple[s]))
    }
    
    tmp_pr_predv[[s]] <- data.frame(pr_coord[[s]] %>% select(lon, lat),tmp_pr_predv[[s]])
    tmp_bg_predv[[s]] <- data.frame(bg_coord[[s]] %>% select(lon, lat),tmp_bg_predv[[s]])
    
    
    eval_res_list[[s]] <- ENMevaluate(occs = tmp_pr_predv[[s]],
                                      bg = tmp_bg_predv[[s]],
                                      tune.args = list(fc = c("L","LQ","LQH","H"),
                                                       rm = c(1,2,4,8,32)),
                                      algorithm = "maxnet",
                                      partitions = "randomkfold",
                                      categoricals = cat,
                                      doClamp = TRUE,
                                      parallel = TRUE,
                                      quiet = TRUE)
    
  }
  
  print("3.models created")
  
  fun_get_AUC_fc <- function(eval_out) {
    opt.aicc <- eval.results(eval_out) %>% filter(delta.AICc == 0)
    data.frame(fc = opt.aicc$fc[1],
               rm = opt.aicc$rm[1],
               AUC = round(opt.aicc$auc.val.avg[1],2))
  }
  
  #output
  df_out$n_pres_pts[c((n*nrow(species_sp)-(nrow(species_sp)-1)):(n*nrow(species_sp)))] <- unlist(lapply(pr_pa, length))
  df_out$n_bg_pts[c((n*nrow(species_sp)-(nrow(species_sp)-1)):(n*nrow(species_sp)))] <- unlist(lapply(bg_pa, length))
  tmp_ind <- which(names(df_out) %in% c("fc", "rm", "AUC"))
  df_out[c((n*nrow(species_sp)-(nrow(species_sp)-1)):(n*nrow(species_sp))),tmp_ind] <- do.call(rbind, lapply(eval_res_list, fun_get_AUC_fc))
  
  
  #2.4. SAC ----
  morans_I_list_p <- list()
  # morans_I_list_pbg <- list()
  for (s in 1:nrow(species_sp)) {
    eval_res <- eval_res_list[[s]]
    res <- eval.results(eval_res)
    opt.aicc <- res %>% filter(delta.AICc == 0)
    mod.seq <- eval.models(eval_res)[[opt.aicc$tune.args[1]]]
    
    #2. make original prediction
    # df_pred_pbg <- rbind(pr_predv[[s]], bg_predv[[s]])
    df_pred_p <- pr_predv[[s]]
    
    # if(all(c("seabed_energy", "seabed_substrate","Month") %in% names(df_pred_pbg))) {
    #   df_pred_pbg <- df_pred_pbg %>%
    #     mutate(sub_char = as.numeric(seabed_substrate)) %>%
    #     mutate(ene_char = as.numeric(seabed_energy)) %>%
    #     mutate(mon_char = as.numeric(Month)) %>%
    #     select(-seabed_substrate, -seabed_energy)
    # }
    if(all(c("seabed_energy", "seabed_substrate","Month") %in% names(df_pred_p))) {
      df_pred_p <- df_pred_p %>%
        mutate(sub_char = as.numeric(seabed_substrate)) %>%
        mutate(ene_char = as.numeric(seabed_energy)) %>%
        mutate(mon_char = as.numeric(Month)) %>%
        select(-seabed_substrate, -seabed_energy)
    }
    # if(all(c("seabed_energy", "seabed_substrate") %in% names(df_pred_pbg))) {
    #   df_pred_pbg <- df_pred_pbg %>%
    #     mutate(sub_char = as.numeric(seabed_substrate)) %>%
    #     mutate(ene_char = as.numeric(seabed_energy)) %>%
    #     select(-seabed_substrate, -seabed_energy)
    # }
    if(all(c("seabed_energy", "seabed_substrate") %in% names(df_pred_p))) {
      df_pred_p <- df_pred_p %>%
        mutate(sub_char = as.numeric(seabed_substrate)) %>%
        mutate(ene_char = as.numeric(seabed_energy)) %>%
        select(-seabed_substrate, -seabed_energy)
    }
    # if(all(c("Month") %in% names(df_pred_pbg))) {
    #   df_pred_pbg <- df_pred_pbg %>%
    #     mutate(mon_char = as.numeric(Month)) %>%
    #     select(-Month)
    # }
    if(all(c("Month") %in% names(df_pred_p))) {
      df_pred_p <- df_pred_p %>%
        mutate(mon_char = as.numeric(Month)) %>%
        select(-Month)
    }
    
    pred_vals_p <- predict(mod.seq, df_pred_p, se.fit=TRUE, type = "cloglog")
    res_p <- c(pr_pa[[s]]) - pred_vals_p
    
    dists_p <- as.matrix(dist(cbind(pr_coord[[s]][,1], pr_coord[[s]][,2])))
    dists.inv_p <- 1/dists_p
    diag(dists.inv_p) <- 0
    dists.inv_p[is.infinite(dists.inv_p)] <- 0 #remove infinite values
    
    # Global Moran's I from ape package
    morans_I_list_p[[s]] <- Moran.I(c(res_p), dists.inv_p, scaled = TRUE, alternative = "greater")
  }
  
  #output
  df_out$SACp[c((n*nrow(species_sp)-(nrow(species_sp)-1)):(n*nrow(species_sp)))] <- round(bind_rows(morans_I_list_p)$observed,3)
  df_out$SACp_pval[c((n*nrow(species_sp)-(nrow(species_sp)-1)):(n*nrow(species_sp)))] <- round(bind_rows(morans_I_list_p)$p.value,5)

  print("4.SAC calculated")
  
  #2.5. response curves ----
  for (s in 1:nrow(species_sp)) {
    #get model
    eval_res <- eval_res_list[[s]]
    opt.aicc <- eval.results(eval_res_list[[s]]) %>% filter(delta.AICc == 0)
    mod <- eval_res@models[[which(names(eval_res@models) == opt.aicc$tune.args[1])]]
    
    plot.new()
    
    #response curves
    png(paste0("SAVE/test_runs_larvae/run_", runX, "/", species_sp$simple[s], fig_out_name, n,".png"), width = 13,
        height = 8, res = 100, units = 'in')
    plot(mod, clamp = T, type = "cloglog")
    mtext(str_to_title(species_sp$simple[s]), side = 3, line = - 2, outer = TRUE)
    dev.off()
  }
  print(paste0("5.loop done - ", n_pres_pts[n], " (run ", runX, ")"))
}

print(df_out)
write.csv(df_out, paste0("SAVE/test_runs_larvae/run_", runX, "/run_results.csv"))
