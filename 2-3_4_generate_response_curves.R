library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)

setwd("DATA")
load("SAVE/NEA_ENMeval_outcomes_adults.Rdata")
load("SAVE/NEA_variable_importance_adults.Rdata")
load("SAVE/NEA_presence_background_adults.Rdata")

pr_predv_ad <- pr_predv
cor_maxnet_ad <- cor_maxnet
eval_res_list_ad <- eval_res_list
rm(bg_coord, bg_pa, bg_predv, pr_coord, pr_pa, pr_predv, cor_maxnet, eval_res_list)

load("SAVE/NEA_ENMeval_outcomes_larvae.Rdata")
load("SAVE/NEA_variable_importance_larvae.Rdata")
load("SAVE/NEA_presence_background_larvae.Rdata")

pr_predv_lv <- pr_predv
rm(bg_coord, bg_pa, bg_predv, pr_coord, pr_pa, pr_predv)


#1. process variable importance ----
var_imp_ad <- cor_maxnet_ad[[1]]
var_imp_lv <- cor_maxnet_sp[[1]]

#convert to percentages
var_imp_ad <- round(var_imp_ad/sum(var_imp_ad),2)
var_imp_lv <- round(var_imp_lv/sum(var_imp_lv),2)

sum(var_imp_ad)
sum(var_imp_lv)

var_imp_ad <- rev(sort(var_imp_ad))
var_imp_lv <- rev(sort(var_imp_lv))


#2. process model outcomes ----

#include only > 5% for plots
var_imp_ad <- var_imp_ad[which(var_imp_ad > 0.05)]
var_imp_lv <- var_imp_lv[which(var_imp_lv > 0.05)]

nms <- union(names(var_imp_ad), names(var_imp_lv))

# create models
eval_res_ad <- eval_res_list_ad[[1]]
res_ad <- eval.results(eval_res_ad)
opt.aicc_ad <- res_ad %>% filter(delta.AICc == 0)
mod_ad <- eval.models(eval_res_ad)[[opt.aicc_ad$tune.args]]

eval_res_lv <- eval_res_list_sp[[1]]
res_lv <- eval.results(eval_res_lv)
opt.aicc_lv <- res_lv %>% filter(delta.AICc == 0)
mod_lv <- eval.models(eval_res_lv)[[opt.aicc_lv$tune.args]]

addline_format <- function(x,...){
  gsub('%','\n',x)
}

name_key <- data.frame(old = c("depth", "sub_char", "SST", "Phyto", "ZooPl", "SSS"),
                       new = c("Depth (m)%", "Seabed substrate", "Sea surface temperature (°C)%", 
                               "Phytoplankton concentration%(mmol C / m³)", "Zooplankton concentration%(g C / m²)", "Sea surface salinity (PSU)%"))
substr_key <- c("Fine mud", "Sandy mud", "Muddy sand", "Sandy mud or Muddy sand ",
                "Sand", "Coarse substrate", "Rock or other hard substrata", 
                "Seabed", "Mixed sediment", "Sediment")

for (n in nms) {
  out_name <- name_key$new[which(name_key$old == n)]
  
  #create a bar plot for seabed substrate
  if(n == "sub_char") {
    dat_lv <- response.plot(mod_lv, n, type = "cloglog", 
                            ylab = "Probability of occurrence",
                            min = min_lv, max = max_lv, plot = F)

    assign(n, ggplot(dat_lv) +
             geom_bar(aes(x = sub_char, y = pred), stat='identity', fill = "#332288") +
             scale_x_discrete(limits = rev(substr_key[which(substr_key %in% dat_lv$sub_char)])) +
             scale_y_continuous(limits = c(0,1)) +
             coord_flip() +
             labs(title = out_name, x = "", y = "Probability of presence") +
             scale_fill_manual(values = c(adult = "#CC6677", larva = "#332288")) + 
             theme_bw() +
             theme(legend.title = element_blank(),
                   plot.title = element_text(size=10, face = "bold", colour = "black", hjust = 0.5),
                   axis.text.y = element_text(size=10, face = "plain", colour = "black"),
                   axis.text.x = element_text(size=10, face = "plain", colour = "black"),
                   axis.title.x = element_text(size=10, face = "bold", colour = "black"),
                   axis.title.y = element_text(size=10, face = "bold", colour = "black"),
                   legend.text = element_text(size=10, face = "bold", colour = "black")))
    next
  }
  
  #the effect of other variables is visualized as a line plot
  
  #calculate x-axis limits
  min_lv <- pr_predv_lv[[1]] %>% dplyr::select(all_of(n)) %>% min
  max_lv <- pr_predv_lv[[1]] %>% dplyr::select(all_of(n)) %>% max
    
  if(n %in% names(var_imp_ad)) {
    min_ad <- pr_predv_ad[[1]] %>% dplyr::select(all_of(n)) %>% min
    max_ad <- pr_predv_ad[[1]] %>% dplyr::select(all_of(n)) %>% max
    
    
    min <- min(min_ad, min_lv)
    max <- max(max_ad, max_lv)
    
    # generate plot data
    dat_ad <- response.plot(mod_ad, n, type = "cloglog", 
                            ylab = "Probability of occurrence",
                            min = min_ad, max = max_ad, plot = F)
    dat_lv <- response.plot(mod_lv, n, type = "cloglog", 
                            ylab = "Probability of occurrence",
                            min = min_lv, max = max_lv, plot = F)
    dat <- rbind(dat_ad %>% mutate(ls = rep("adult")), 
                 dat_lv %>% mutate(ls = rep("larva")))
    
    # plot using ggplot2
    assign(n, ggplot(dat) +
      geom_line(aes_string(x = n, y = "pred", colour = "ls"), linewidth = 1) + 
      labs(x = addline_format(out_name), y = "Probability of presence") +
      scale_colour_manual(values = c(adult = "#CC6677", larva = "#332288"), name = "Life stage") + 
      scale_x_continuous(limits = c(min, max)) +
      scale_y_continuous(limits = c(0,1)) +
      theme_bw() +
      theme(legend.title = element_text(size=10, face = "bold", colour = "black"),
            legend.text = element_text(size=10, face = "plain", colour = "black"), 
            axis.text.y = element_text(size=10, face = "plain", colour = "black"),
            axis.text.x = element_text(size=10, face = "plain", colour = "black"),
            axis.title.x = element_text(size=10, face = "bold", colour = "black"),
            axis.title.y = element_text(size=10, face = "bold", colour = "black")))
  
  }
  
  #since larvae have more variables that have var importance > 5 %, some plots need to be generated only for larvae
  else if(!(n %in% names(var_imp_ad))) {
    dat_lv <- response.plot(mod_lv, n, type = "cloglog", 
                            ylab = "Probability of occurrence",
                            min = min_lv, max = max_lv, plot = F)
    
    assign(n, ggplot() +
             geom_line(data = dat_lv, aes_string(x = n, y = "pred"), colour = "#332288", linewidth = 1) + 
             labs(x = addline_format(out_name), y = "Probability of presence") +
             scale_x_continuous(limits = c(min_lv, max_lv)) +
             scale_y_continuous(limits = c(0,1)) +
             theme_bw() +
             theme(legend.title = element_blank(),
                   axis.text.y = element_text(size=10, face = "plain", colour = "black"),
                   axis.text.x = element_text(size=10, face = "plain", colour = "black"),
                   axis.title.x = element_text(size=10, face = "bold", colour = "black"),
                   axis.title.y = element_text(size=10, face = "bold", colour = "black")))
    }
}

# compile all single plots in one figure

# generate common legend:
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(depth + theme(legend.background = element_rect(linetype = 1, linewidth = 0.1, colour = 1)))


a <- grid.arrange(depth + theme(legend.position = "none"), 
                  SST + labs(y = "") + theme(axis.text.y=element_blank(), legend.position = "none"), 
                  SSS + labs(y = "") + theme(axis.text.y=element_blank(), legend.position="none"), 
                  ZooPl + theme(legend.position = "none"), 
                  Phyto + labs(y = "") + theme(axis.text.y=element_blank(), legend.position = "none"), 
                  sub_char,
                  mylegend ,
                  ncol = 6,
                  layout_matrix = rbind(c(1,1,2,2,3,3),
                                        c(4,4,5,5,7,7),
                                        c(6,6,6,6,7,7)))

# save plot
ggsave("3.MODEL_OUTPUT/response_curves.png", a, width = 10.4, height = 8.66)




