# Load data ---------------------------------------------------------------

# Loop through all the csv files and load them into R
# !!! CHANGE THIS PATH !!!
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018'

file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'breeder_stats_bydiam_2census', 'breederscore_bin_bydiam_2census', 'breeder_stats_bylight_2census', 'breederscore_bin_bylight_2census', 'fastslow_stats_bydiam_2census', 'fastslowscore_bin_bydiam_2census', 'fastslow_stats_bylight_2census', 'fastslowscore_bin_bylight_2census')

for (i in file_names) {
  assign(i, read.csv(file.path(fpdata, paste0(i,'.csv')), stringsAsFactors = FALSE))
}

# Plot data ---------------------------------------------------------------

library(cowplot)
library(dplyr)

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.

guild_colors <- RColorBrewer::brewer.pal(5, 'Set1')
fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('fast','long-lived pioneer', 'slow', 'short-lived breeder', 'intermediate')

# Figure 3a
fig_3a <- indivproductionbin_5census %>%
  filter(fg %in% fg_names) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

# Figure 3b
fig_3b <- densitybin_5census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

# Figure 3c
fig_3c <- totalproductionbin_5census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Total production (kg ha'^-1,' y'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

# Figure 4a
fig_4a <- crownareabin_2census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Total crown area (m'^2, ' ha'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

# Figure 4b
fig_4b <- lightreceivedbin_2census %>%
  filter(fg %in% fg_names) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = expression(paste('Total light received (W ha'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

# Figure 4c
fig_4c <- indivprodperareabin_2census %>%
  filter(fg %in% fg_names) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = expression(paste('Production per unit crown area (kg y'^-1, ' m'^-2,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels) +
  panel_border(colour = 'black')

# Make a figure to show the FGs
ggplot(fgbci, aes(x = X1, y = X2, color = factor(fg))) +
  geom_point() +
  labs(x = 'X1 slow to fast', y = 'X2 pioneers to breeders') +
  scale_color_manual(values = guild_colors)

###
# Added 23 Jan: Ratio figures (fig 5)
# Note that the continuous figures, as they don't use groupings, use all the data.

# Breeder to pioneer by light
fig_blightprod <- breeder_stats_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Breeder to pioneer production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_blightdens <- breeder_stats_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Breeder to pioneer density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_blightscore <- breederscore_bin_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_continuous(name = 'Low = long-lived pioneer, high = short-lived breeder') +
  panel_border(colour = 'black')

# Breeder to pioneer by diameter
fig_bdiamprod <- breeder_stats_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Breeder to pioneer production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_bdiamdens <- breeder_stats_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Breeder to pioneer density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_bdiamscore <- breederscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_continuous(name = 'Low = long-lived pioneer, high = short-lived breeder') +
  panel_border(colour = 'black')


####

# Fast to slow by light
fig_flightprod <- fastslow_stats_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Fast to slow production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_flightdens <- fastslow_stats_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Fast to slow density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_flightscore <- fastslowscore_bin_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_continuous(name = 'Low = slow, high = fast') +
  panel_border(colour = 'black')

# Fast to slow by diameter
fig_fdiamprod <- fastslow_stats_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Fast to slow production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_fdiamdens <- fastslow_stats_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_log10(name = 'Fast to slow density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_fdiamscore <- fastslowscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)') + 
  scale_y_continuous(name = 'Low = slow, high = fast') +
  panel_border(colour = 'black')
