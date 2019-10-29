# Load data ---------------------------------------------------------------

# Loop through all the csv files and load them into R
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_june2018_alternativecluster'
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'breeder_stats_bydiam_2census', 'breederscore_bin_bydiam_2census', 'breeder_stats_bylight_2census', 'breederscore_bin_bylight_2census', 'fastslow_stats_bydiam_2census', 'fastslowscore_bin_bydiam_2census', 'fastslow_stats_bylight_2census', 'fastslowscore_bin_bylight_2census')

for (i in file_names) {
  assign(i, read.csv(file.path(fpdata, paste0(i,'.csv')), stringsAsFactors = FALSE))
}

# Plot data ---------------------------------------------------------------

library(cowplot)
library(dplyr)

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.
global_diam_xlimits <- c(1, 316) # Maximum and minimum for x axis if it's diameter.

guild_colors <- RColorBrewer::brewer.pal(5, 'Set1')
fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('fast','long-lived pioneer', 'slow', 'short-lived breeder', 'intermediate')

# Figure 3a
fig_3a <- indivproductionbin_5census %>%
  filter(fg %in% fg_names, !is.na(mean), mean > 0) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group') +
  panel_border(colour = 'black')

# Figure 3b
fig_3b <- densitybin_5census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 3c
fig_3c <- totalproductionbin_5census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Total production (kg ha'^-1,' y'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 4a
fig_4a <- crownareabin_2census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Total crown area (m'^2, ' ha'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 4b
fig_4b <- lightreceivedbin_2census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Total light received (W ha'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 4c
fig_4c <- indivprodperareabin_2census %>%
  filter(fg %in% fg_names, !is.na(mean), mean > 0) %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = expression(paste('Production per unit crown area (kg y'^-1, ' m'^-2,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group') +
  panel_border(colour = 'black')

pdf('C:/Users/Q/google_drive/ForestLight/figs/figures_june2018_newcluster/figs3and4.pdf', height = 6, width = 7.5)
print(fig_3a + ggtitle('3A'))
print(fig_3b + ggtitle('3B'))
print(fig_3c + ggtitle('3C'))
print(fig_4a + ggtitle('4A'))
print(fig_4b + ggtitle('4B'))
print(fig_4c + ggtitle('4C'))
dev.off()

# Plot functional groups
ggplot(fgbci, aes(x = PC_slow_to_fast, y = PC_breeder_to_pioneer, color = factor(fg5))) +
  geom_point() +
  labs(x = 'X1 slow to fast', y = 'X2 breeders to pioneers') +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group')
ggsave('C:/Users/Q/google_drive/ForestLight/figs/figures_june2018_newcluster/fg5plot (newest groups).pdf')
###
# Added 23 Jan: Ratio figures (fig 5)
# Note that the continuous figures, as they don't use groupings, use all the data.

# Breeder to pioneer by light
fig_blightprod <- breeder_stats_bylight_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Breeder to pioneer production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_blightdens <- breeder_stats_bylight_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
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
  scale_y_continuous(name = 'Low = short-lived breeder, high = long-lived pioneer') +
  panel_border(colour = 'black')

# Breeder to pioneer by diameter
fig_bdiamprod <- breeder_stats_bydiam_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = 'Breeder to pioneer production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_bdiamdens <- breeder_stats_bydiam_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = 'Breeder to pioneer density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_bdiamscore <- breederscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_continuous(name = 'Low = short-lived breeder, high = long-lived pioneer') +
  panel_border(colour = 'black')


####

# Fast to slow by light
fig_flightprod <- fastslow_stats_bylight_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Fast to slow production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_flightdens <- fastslow_stats_bylight_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
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
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = 'Fast to slow production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_fdiamdens <- fastslow_stats_bydiam_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = 'Fast to slow density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_fdiamscore <- fastslowscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_continuous(name = 'Low = slow, high = fast') +
  panel_border(colour = 'black')

pdf('C:/Users/Q/google_drive/ForestLight/figs/figures_june2018_newcluster/ratio_figures.pdf', height = 6, width = 6)
print(fig_blightdens)
print(fig_blightprod)
print(fig_blightscore)
print(fig_bdiamdens)
print(fig_bdiamprod)
print(fig_bdiamscore)
print(fig_flightdens)
print(fig_flightprod)
print(fig_flightscore)
print(fig_fdiamdens)
print(fig_fdiamprod)
print(fig_fdiamscore)
dev.off()