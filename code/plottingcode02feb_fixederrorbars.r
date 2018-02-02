# PLOTTING CODE 02 FEB
# NOTE THAT DIFFERENT DATA IS LOADED DEPENDING ON WHETHER YOU DON'T WANT NO SHRUB
# EDITED VERSION: 
# Edit 1: Error bars have horizontal caps with fixed width
# Edit 2: At the bottom of this script, alternative versions of the quantile plots are made with 95% CI of the mean instead

include_shrubs <- TRUE # *** change this to FALSE to exclude shrubs.

# Load data ---------------------------------------------------------------

# Loop through all the csv files and load them into R
fpdata <- 'C:/Users/Q/google_drive/ForestLight/data/data_22jan2018' # *** change this path to correct path
file_names <- c('densitybin_5census', 'indivproductionbin_5census', 'totalproductionbin_5census', 'crownareabin_2census', 'lightreceivedbin_2census', 'indivprodperareabin_2census', 'breeder_stats_bydiam_2census', 'breederscore_bin_bydiam_2census', 'breeder_stats_bylight_2census', 'breederscore_bin_bylight_2census', 'fastslow_stats_bydiam_2census', 'fastslowscore_bin_bydiam_2census', 'fastslow_stats_bylight_2census', 'fastslowscore_bin_bylight_2census')

fgbci <- read.table('C:/Users/Q/google_drive/ForestLight/data/Ruger/functional_groups_BCI.txt', stringsAsFactors = FALSE)

if (include_shrubs) {
  for (i in file_names) {
    assign(i, read.csv(file.path(fpdata, paste0(i,'.csv')), stringsAsFactors = FALSE))
  }
} else {
  for (i in file_names) {
    assign(i, read.csv(file.path(fpdata, paste0('noshrub_', i,'.csv')), stringsAsFactors = FALSE))
    fgbci <- subset(fgbci, !grform %in% 'S')
  }
}

# Plot data ---------------------------------------------------------------

library(cowplot)
library(dplyr)

# Set options for error bar widths and dodge amounts for all the plots
error_bar_width <- 0.03
dodge_width <- 0.03

area_core <- 42.84 # Area of the plot without the edge strip and without the "young" forest area.
global_diam_xlimits <- c(1, 316) # Maximum and minimum for x axis if it's diameter.

guild_colors <- RColorBrewer::brewer.pal(5, 'Set1')
fg_names <- paste('fg', 1:5, sep = '')
fg_labels <- c('fast','long-lived pioneer', 'slow', 'short-lived breeder', 'intermediate')

p_dodge <- position_dodge(width = dodge_width)

# Figure 3a
fig_3a <- indivproductionbin_5census %>%
  filter(fg %in% fg_names, !is.na(mean), mean > 0) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_errorbar(aes(width = width), position=p_dodge) +
  geom_point(position=p_dodge) +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group') +
  panel_border(colour = 'black')

# Figure 3b
fig_3b <- densitybin_5census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_errorbar(aes(width = width), position=p_dodge) +
  geom_point(position=p_dodge) +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Density (individuals ha'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 3c
fig_3c <- totalproductionbin_5census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_errorbar(aes(width = width), position=p_dodge) +
  geom_point(position=p_dodge) +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Total production (kg ha'^-1,' y'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 4a
fig_4a <- crownareabin_2census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_errorbar(aes(width = width), position=p_dodge) +
  geom_point(position=p_dodge) +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Total crown area (m'^2, ' ha'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 4b
fig_4b <- lightreceivedbin_2census %>%
  filter(fg %in% c('all', fg_names), !is.na(bin_yvalue), bin_yvalue > 0) %>%
  mutate(bin_ymin = ifelse(bin_ymin == 0, bin_yvalue, bin_ymin)) %>%
  mutate(bin_yvalue = bin_yvalue/area_core, bin_ymin = bin_ymin/area_core, bin_ymax = bin_ymax/area_core) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = bin_yvalue, ymin = bin_ymin, ymax = bin_ymax, group = fg, color = fg)) +
  geom_errorbar(aes(width = width), position=p_dodge) +
  geom_point(position=p_dodge) +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Total light received (W ha'^-1,')'))) +
  scale_color_manual(values = c('black', guild_colors), labels = c('all', fg_labels), name = 'functional group') +
  panel_border(colour = 'black')

# Figure 4c
fig_4c <- indivprodperareabin_2census %>%
  filter(fg %in% fg_names, !is.na(mean), mean > 0) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975, group = fg, color = fg)) +
  geom_errorbar(aes(width = width)) +
  geom_point() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = expression(paste('Production per unit crown area (kg y'^-1, ' m'^-2,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group') +
  panel_border(colour = 'black')

# Plot functional groups
ggplot(fgbci, aes(x = X1, y = X2, color = factor(fg))) +
  geom_point() +
  labs(x = 'X1 slow to fast', y = 'X2 pioneers to breeders') +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group')


# Plot ratios -------------------------------------------------------------

# Note that the continuous figures, as they don't use groupings, use all the data.

# Breeder to pioneer by light
fig_blightprod <- breeder_stats_bylight_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_errorbar(width = error_bar_width) +
  geom_point() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Breeder to pioneer production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_blightdens <- breeder_stats_bylight_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_errorbar(width = error_bar_width) +
  geom_point() +
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
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_errorbar(width = error_bar_width) +
  geom_point() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = 'Breeder to pioneer production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_bdiamdens <- breeder_stats_bydiam_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_errorbar(width = error_bar_width) +
  geom_point() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = 'Breeder to pioneer density ratio', breaks = c(0.1,1,10,100)) +
  panel_border(colour = 'black')

fig_bdiamscore <- breederscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = q025, ymax = q975)) +
  geom_segment(aes(xend = bin_midpoint, y = q25, yend = q75), size = 2, color = 'gray50') +
  geom_pointrange() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_continuous(name = 'Low = long-lived pioneer, high = short-lived breeder') +
  panel_border(colour = 'black')


####

# Fast to slow by light
fig_flightprod <- fastslow_stats_bylight_2census %>%
  filter(production_ratio_mean > 0) %>%
  mutate(production_ratio_min = ifelse(production_ratio_min == 0, production_ratio_mean, production_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = production_ratio_mean, ymin = production_ratio_min, ymax = production_ratio_max)) +
  geom_errorbar(width = error_bar_width) +
  geom_point() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = 'Fast to slow production ratio', breaks= 10^(-1:3)) +
  panel_border(colour = 'black')

fig_flightdens <- fastslow_stats_bylight_2census %>%
  filter(density_ratio_mean > 0) %>%
  mutate(density_ratio_min = ifelse(density_ratio_min == 0, density_ratio_mean, density_ratio_min)) %>%
  ggplot(aes(x = bin_midpoint, y = density_ratio_mean, ymin = density_ratio_min, ymax = density_ratio_max)) +
  geom_errorbar(width = error_bar_width) +
  geom_point() +
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
  geom_errorbar(width = error_bar_width) +
  geom_point() +
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

#############################################################################################
#############################################################################################
# Everything repeated with means and confidence intervals instead of quantiles
# (or at least the ones that use quantiles above)
#############################################################################################
#############################################################################################


# Figure 3a
fig_3a_CI <- indivproductionbin_5census %>%
  filter(fg %in% fg_names, !is.na(mean), mean > 0) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, group = fg, color = fg)) +
  geom_errorbar(aes(width = width), position=p_dodge) +
  geom_point(position=p_dodge) +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_log10(name = expression(paste('Production (kg y'^-1,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group') +
  panel_border(colour = 'black')

# Figure 4c
fig_4c_CI <- indivprodperareabin_2census %>%
  filter(fg %in% fg_names, !is.na(mean), mean > 0) %>%
  group_by(bin_midpoint) %>% mutate(width = error_bar_width * n()) %>% ungroup %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max, group = fg, color = fg)) +
  geom_errorbar(aes(width = width)) +
  geom_point() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_log10(name = expression(paste('Production per unit crown area (kg y'^-1, ' m'^-2,')'))) +
  scale_color_manual(values = guild_colors, labels = fg_labels, name = 'functional group') +
  panel_border(colour = 'black')

# Breeder to pioneer by light
fig_blightscore_CI <- breederscore_bin_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max)) +
  geom_errorbar(width = error_bar_width) +
  geom_point() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_continuous(name = 'Low = long-lived pioneer, high = short-lived breeder') +
  panel_border(colour = 'black')

# Breeder to pioneer by diameter
fig_bdiamscore_CI <- breederscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max)) +
  geom_errorbar(width = error_bar_width) +
  geom_point() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_continuous(name = 'Low = long-lived pioneer, high = short-lived breeder') +
  panel_border(colour = 'black')


####

# Fast to slow by light
fig_flightscore_CI <- fastslowscore_bin_bylight_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max)) +
  geom_errorbar(width = error_bar_width) +
  geom_point() +
  scale_x_log10(name = expression(paste('Light received per unit crown area (W m'^-2,')'))) + 
  scale_y_continuous(name = 'Low = slow, high = fast') +
  panel_border(colour = 'black')

# Fast to slow by diameter
fig_fdiamscore_CI <- fastslowscore_bin_bydiam_2census %>%
  ggplot(aes(x = bin_midpoint, y = mean, ymin = ci_min, ymax = ci_max)) +
  geom_errorbar(width = error_bar_width) +
  geom_point() +
  scale_x_log10(name = 'Diameter (cm)', limits = global_diam_xlimits) + 
  scale_y_continuous(name = 'Low = slow, high = fast') +
  panel_border(colour = 'black')
