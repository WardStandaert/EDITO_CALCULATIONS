library(zoo)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)

setwd("DATA")

df_HSI_raw <- read.csv("3.MODEL_OUTPUT/LARVAE/HSI_variation_BPNS.csv")
df_NAO_raw <- read.delim("1.DOWNLOAD/NAO_indices/nao_station_monthly.txt", header = T, skip = 1, sep = "")

#define seasons
seasons <- data.frame(month = month.abb,
                      season = rep(c("winter","spring","summer","autumn"), each = 3))


#1. prepare autumn NAO 5-year average ----

#prepare NAO data frame
df_NAO <- df_NAO_raw %>% 
  filter(as.numeric(row.names(.)) >= 1998 & as.numeric(row.names(.)) <= 2022) %>%   #keep only data from 1998 - 2022
  mutate(year = row.names(.)) %>%
  gather(key = month, value = NAO, -year) %>%                                       #convert format
  left_join(seasons, by = "month") %>%                                              #add seasons
  group_by(year, season) %>% 
  dplyr::summarize(NAO = mean(NAO)) %>%                                             #average NAO per year & season
  mutate(year = as.numeric(year))

#calculate 5 year average for the autumn NAO
df_NAO_5y_av <- df_NAO %>% 
  filter(season == "autumn") %>%
  select(-season) %>%
  zoo::rollmean(5) %>%
  as.data.frame()


#2. prepare winter HSI 5-year average ----

#calculate 5 year average for larval HSI (best months for BPNS: December & January)
df_HSI <- df_HSI_raw %>%
  select(-X) %>%
  gather(key = year, value = HSI, -month) %>%                   #convert format
  filter(month %in% c("Dec","Jan")) %>%                         #keep only December & January
  mutate(year2 = as.numeric(str_remove(year,"X"))) %>%          #convert year into numeric values
  mutate(year = ifelse(month == "Jan", year2 - 1, year2)) %>%   #align combine January of next year with December of previous year (now year refers to starting year of the winter)
  group_by(year) %>% 
  dplyr::summarize(HSI = mean(HSI)) %>%                         #average HSI per year
  filter(year >= 2000)                                          #remove 1999 (no values for December 1999)

#calculate 5 year average
df_HSI_5y_av <- df_HSI %>%
  rollmean(5) %>%
  as.data.frame()


#3. combine, correlate and plot ----
df_plot <- left_join(df_HSI_5y_av, df_NAO_5y_av, by = "year")

cor(df_plot$NAO, df_plot$HSI)
# -0.7211571

#define limits for primary and secondary y-axis
ylim_prim <- c(0.4, 0.7)
ylim_sec  <- c(1.2, -0.6)

b <- diff(ylim_prim)/diff(ylim_sec)
a <- ylim_prim[1] - b*ylim_sec[1]

p1 <- ggplot(df_plot) +
  geom_line(aes(x = year, y = HSI), linewidth = .8) +
  geom_line(aes(x = year, y = a + NAO*b), color = "#CC6677", linewidth = .8) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by = 0.1), 
                     limits = c(0.4, 0.7),
                     name = "Winter HSI",
                     sec.axis = sec_axis(~ (. - a)/b, name = "Autumn NAO", breaks = seq(-0.6, 1.2, 0.3))) +
  scale_x_continuous(breaks = seq(2000,2020,2)) +
  theme_bw() +
  theme(panel.grid.minor =  element_blank(),
        panel.grid.major.x =  element_blank(),
        axis.text.y = element_text(size=9, face = "plain", colour = "black"),
        axis.text.x = element_text(size=9, face = "plain", colour = "black"),
        axis.title.x = element_text(size=10, face = "bold", colour = "black"),
        axis.title.y = element_text(size=10, face = "bold", colour = "black"),
        axis.line.y.right = element_line(color = "#CC6677"), 
        axis.ticks.y.right = element_line(color = "#CC6677"),
        axis.text.y.right = element_text(color = "#CC6677"), 
        axis.title.y.right = element_text(color = "#CC6677"))

p2 <- ggplot(df_plot) +
  geom_point(aes(x = NAO, y = HSI)) +
  scale_y_continuous(breaks = seq(0.4,0.7,0.1), limits = c(0.4,0.7))+
  geom_vline(xintercept = 0, linetype = 3) +
  xlab("Autumn NAO") + 
  ylab("Winter HSI") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size=9, face = "plain", colour = "black"),
        axis.text.x = element_text(size=9, face = "plain", colour = "black"),
        axis.title.x = element_text(size=10, face = "bold", colour = "black"),
        axis.title.y = element_text(size=10, face = "bold", colour = "black"))

#combine plots and save
ggarrange(p1, p2, labels = c("a.","b."), font.label = list(size = 10, color = "black", face = "bold"))

ggsave("3.MODEL_OUTPUT/LARVAE/HSI_vs_NAO.png", width = 9, height = 3)
