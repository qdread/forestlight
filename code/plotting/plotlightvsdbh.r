# Light Vs Size Plot

gdrive_path <- '~/google_drive'
github_path <- '~/Documents/GitHub'

load(file.path(gdrive_path, 'ForestLight/data/rawdataobj_alternativecluster.r'))
source(file.path(github_path, 'forestlight/code/allfunctions27july.r'))

library(dplyr)
library(ggplot2)

# Do some bins ------------------------------------------------------------


alltree_light_95 <- alltree_light_95 %>%
  mutate(fg = if_else(!is.na(fg), paste0('fg',fg), as.character(NA)))

dbhbin1995 <- with(alltree_light_95, logbin(x = dbh_corr, n = 20))

lightbydbhfakebin_fg <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(fg, dbh_bin) %>%
  do(quantile(.$light_received/.$crownarea, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))

lightbydbhfakebin_all <- alltree_light_95 %>%
  mutate(dbh_bin = cut(dbh_corr, breaks = c(dbhbin1995$bin_min[1], dbhbin1995$bin_max), labels = dbhbin1995$bin_midpoint, include.lowest = TRUE)) %>%
  group_by(dbh_bin) %>%
  do(quantile(.$light_received/.$crownarea, c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% 
       as.data.frame %>%
       setNames(c('q025','q25','q50','q75','q975')))
  
lightbydbhfakebin_fg <- data.frame(fg = 'all', lightbydbhfakebin_all, stringsAsFactors = FALSE) %>%
  rbind(as.data.frame(lightbydbhfakebin_fg)) %>%
  mutate(dbh_bin = as.numeric(as.character(dbh_bin)))


# Plot: raw data ----------------------------------------------------------

exl <- expression(paste('Light received per crown area (W m'^-2, ')', sep = ''))
exd <- 'Diameter (cm)'

# Each group
ggplot() +
  geom_point(alpha = 0.05, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received/crownarea), color = 'skyblue') +
  geom_pointrange(data = lightbydbhfakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_bw() 

# All together
ggplot() +
  geom_point(alpha = 0.05, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownarea), color = 'skyblue') +
  geom_pointrange(data = lightbydbhfakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_bw() 


# Plot: hexagon plot ------------------------------------------------------

alpha_value <- 0.6
hexfill <- scale_fill_gradient(low = 'skyblue', high = 'navy', trans = 'log', breaks = c(1,3,10,30,100,300))

# Each group
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95 %>% filter(!is.na(fg)), aes(x = dbh_corr, y = light_received/crownarea)) +
  geom_pointrange(data = lightbydbhfakebin_fg %>% filter(!fg %in% 'all', !is.na(fg)), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  facet_wrap(~ fg) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_bw() +
  hexfill +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

# All together
ggplot() +
  geom_hex(alpha = alpha_value, data = alltree_light_95, aes(x = dbh_corr, y = light_received/crownarea)) +
  geom_pointrange(data = lightbydbhfakebin_fg %>% filter(fg %in% 'all'), aes(x = dbh_bin, y = q50, ymin = q25, ymax = q75)) +
  scale_x_log10(name = exd) +
  scale_y_log10(name = exl) +
  theme_bw() +
  hexfill +
  guides(fill = guide_legend(override.aes = list(alpha = alpha_value)))

