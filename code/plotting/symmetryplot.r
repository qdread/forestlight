# Symmetry plot
# Compare total growth fitted slope (evaluated halfway between the high and low cutoffs) with total light fitted slope at the same point. 
# QDR / Forestlight / 20 June 2019

library(tidyverse)

gdrive_path <- '~/google_drive/ForestLight'

# Load parameter df to find the point at which to evaluate the fitted slopes, then load the fitted slopes

params <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/newpiecewise_paramci_by_fg.csv'), stringsAsFactors = FALSE)

xvalues <- params %>% 
  filter(variable == 'density', model == 3) %>%
  group_by(fg) %>%
  summarize(mid_cutoff = (mean[parameter == 'tau_high'] + mean[parameter == 'tau_low'])/2)

# Load fitted slopes of total growth and total light.
growth_slopes <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/newpiecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)
light_slopes <- read.csv(file.path(gdrive_path, 'data/data_piecewisefits/totallightscaling/light_piecewise_fitted_slopes_by_fg.csv'), stringsAsFactors = FALSE)

growth_slopes_atmiddle <- growth_slopes %>% 
  filter(variable == 'total_production', dens_model == 3, prod_model == 2) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

light_slopes_atmiddle <- light_slopes %>% 
  filter(variable == 'total_incoming_light', dens_model == 3, prod_model == 2) %>%
  left_join(xvalues) %>%
  mutate(diff = abs(dbh - mid_cutoff)) %>%
  group_by(fg) %>%
  filter(diff == min(diff))

fg_full_names <- c('fast', 'LL pioneer', 'slow', 'SL breeder', 'medium', 'all trees', 'unclassified')
fgs <- c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified')

allslopes <- rbind(growth_slopes_atmiddle, light_slopes_atmiddle) %>%
  ungroup %>%
  mutate(fg = factor(fg, levels = fgs, labels = fg_full_names))
  

ggplot(allslopes %>% filter(!fg %in% 'unclassified'), aes(x = fg, y = q50, ymin = q025, ymax = q975, color = variable)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_errorbar(position = position_dodge(width = 0.05), size = 1, width = 0.2) +
  geom_point(position = position_dodge(width = 0.05), size = 2) +
  labs(x = 'functional group', y = 'scaling slope for midsize trees') +
  scale_color_manual(values = c('slateblue', 'goldenrod'), labels = c('total light', 'total growth')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom')
