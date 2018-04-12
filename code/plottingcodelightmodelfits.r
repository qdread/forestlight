# Plotting code for production(growth) per area versus incoming light per area
# Plots the functional fits over the raw data
# Also plots the functional fits over the medians and quantiles, and the means and confidence intervals of each bin.

year_to_plot <- 1995 ### CHANGE THIS IF YOU WANT TO PLOT 1990

# Load data ---------------------------------------------------------------

fp <- 'C:/Users/Q/google_drive/ForestLight/data/data_forplotting_light_12apr2018'

obs_light_binned <- read.csv(file.path(fp, 'obs_light_binned.csv'), stringsAsFactors = FALSE)
obs_light_raw <- read.csv(file.path(fp, 'obs_light_raw.csv'), stringsAsFactors = FALSE)
pred_light <- read.csv(file.path(fp, 'pred_light.csv'), stringsAsFactors = FALSE)
param_ci <- read.csv(file.path(fp, 'lightbyarea_paramci_by_fg.csv'), stringsAsFactors = FALSE)

library(dplyr)
library(cowplot)

# Get rid of the predicted points that are outside the limits of the observed data for each FG
obs_limits <- obs_light_binned %>%
  group_by(fg, year) %>%
  summarize(min_obs = min(bin_midpoint), max_obs = max(bin_midpoint))

pred_light <- pred_light %>%
  left_join(obs_limits) %>%
  filter(light_area >= min_obs & light_area <= max_obs)

pred_light_5groups <- pred_light %>% filter(!fg %in% c('alltree','unclassified'))

# Do some additional computation to correct the error bar width for the number of groups in each bin
obs_light_binned <- obs_light_binned %>%
  group_by(bin_midpoint, year) %>% mutate(width = sum(c('fg1','fg2','fg3','fg4','fg5') %in% fg)) %>% ungroup

# Create plots ------------------------------------------------------------

### ------------------------------------------------------------------------------------------
### SET THESE OPTIONS

# Names of functional groups to display
fg_display <- c(fg1 = 'fast', fg2 = 'slow', fg3 = 'pioneer', fg4 = 'breeder', fg5 = 'middle')

# Axis titles
title_x <- expression(paste('Incoming light per area (W m'^-2,')',sep=''))
title_y <- expression(paste('Growth per area (kg y'^-1, ' m'^-2,')', sep=''))

# Colors
# these are shit colors but I just put them in as a placeholder

fg_colors <- c(fg1 = 'red', fg2 = 'blue', fg3 = 'green', fg4 = 'purple', fg5 = 'orange', unclassified = 'brown', alltree = 'black')

### -----------------------------------------------------------------------------------------

# 1. Plot of maximum slope by functional group

# Remove all tree and unclassified groups

ggplot(param_ci %>% filter(year == year_to_plot, parameter %in% 'max_slope', !fg %in% c('alltree','unclassified')),
       aes(x = fg, y = q50, ymin = q025, ymax = q975)) +
  geom_errorbar(width = 0.1) + geom_point() +
  labs(y = 'maximal slope') +
  scale_x_discrete(name = 'functional group', labels = fg_display) +
  panel_border(colour = 'black')

# 2. Plot with different panels for each functional group, and raw data

# I attempted to set an alpha scale so that the amount of transparency is roughly the same but the numbers may need to be tweaked

p_raw_panels <- ggplot(obs_light_raw %>% filter(year == year_to_plot, !fg %in% 'unclassified')) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_display)) +
  geom_point(aes(x = light_area, y = production_area, alpha = fg)) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975), alpha = 0.5, fill = 'red') +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50), color = 'red') +
  scale_alpha_manual(values = c(0.1, 0.08, 0.02, 0.02, 0.004)) +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA),
        legend.position = 'none')

# 3. Plot with different panels for each functional group, and quantiles

p_median_panels <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_display)) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975), alpha = 0.5, fill = 'red') +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50), color = 'red') +
  geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q25, yend = q75), size = 0.75) +
  geom_segment(aes(x = bin_midpoint, xend = bin_midpoint, y = q025, yend = q975)) +
  geom_point(aes(x = bin_midpoint, y = median)) +
  
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA))

# 4. Plot with different panels for each functional group, and means

p_mean_panels <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  facet_wrap(~ fg, labeller = labeller(fg = fg_display)) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975), alpha = 0.5, fill = 'red') +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50), color = 'red') +
  geom_errorbar(aes(x = bin_midpoint, ymin = ci_min, ymax = ci_max)) +
  geom_point(aes(x = bin_midpoint, y = mean)) +
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA))

# 5. Plot with all functional groups on the same panel, and quantiles

dodge_width <- 0.03

p_median_1panel <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50, color = fg)) +
  geom_errorbar(aes(x = bin_midpoint, ymin = q25, ymax = q75, group = fg, color = fg), size = 0.75, width = 0, position = position_dodge(width = dodge_width)) +
  geom_errorbar(aes(x = bin_midpoint, ymin = q025, ymax = q975, group = fg, color = fg), width = 0, position = position_dodge(width = dodge_width)) +
  geom_point(aes(x = bin_midpoint, y = median, group = fg, color = fg), position = position_dodge(width = dodge_width)) +
  
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y) +
  scale_color_manual(name = 'Functional group', values = fg_colors, labels = fg_display) +
  scale_fill_manual(values = fg_colors, labels = fg_display, guide = FALSE) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        legend.position = c(0.2, 0.8))

# 6. Plot with all functional groups on the same panel, and means

dodge_width <- 0.03
error_bar_width <- 0.03

p_mean_1panel <- ggplot(obs_light_binned %>% filter(year == year_to_plot, !fg %in% c('alltree', 'unclassified'))) +
  geom_ribbon(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, ymin = q025, ymax = q975, fill = fg), alpha = 0.5) +
  geom_line(data = pred_light_5groups %>% filter(year == year_to_plot), aes(x = light_area, y = q50, color = fg)) +
  geom_errorbar(aes(x = bin_midpoint, ymin = ci_min, ymax = ci_max, group = fg, color = fg, width = error_bar_width * width), position = position_dodge(width = dodge_width)) +
  geom_point(aes(x = bin_midpoint, y = mean, group = fg, color = fg), position = position_dodge(width = dodge_width)) +
  
  scale_x_log10(name = title_x) + 
  scale_y_log10(name = title_y) +
  scale_color_manual(name = 'Functional group', values = fg_colors, labels = fg_display) +
  scale_fill_manual(values = fg_colors, labels = fg_display, guide = FALSE) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA),
        legend.position = c(0.2, 0.8))


